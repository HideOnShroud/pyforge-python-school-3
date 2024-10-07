from airflow import DAG
from airflow.operators.python import PythonOperator
from airflow.utils.dates import days_ago
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine
import pandas as pd
from dotenv import load_dotenv
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
import boto3
from datetime import datetime
import os


DATABASE_URL = "postgresql://admin:password@postgres:5432/drug_db"
engine = create_engine(DATABASE_URL)
SessionLocal = sessionmaker(bind=engine)

load_dotenv("/.env")

# Load S3 configurations
S3_BUCKET = os.getenv("S3_BUCKET")
S3_ACCESS_KEY = os.getenv("S3_ACCESS_KEY")
S3_SECRET_KEY = os.getenv("S3_SECRET_KEY")
S3_ENDPOINT = os.getenv("S3_ENDPOINT")

# Initialize Minio client or boto3 S3
s3_client = boto3.client(
    's3',
    endpoint_url=S3_ENDPOINT,
    aws_access_key_id=S3_ACCESS_KEY,
    aws_secret_access_key=S3_SECRET_KEY
)

# Default arguments for the DAG
default_args = {
    'owner': 'airflow',
    'depends_on_past': False,
    'start_date': days_ago(1),
    'email_on_failure': False,
    'email_on_retry': False
}

# Define the DAG
dag = DAG(
    'molecule_etl_dag',
    default_args=default_args,
    description='Extract, transform and load molecules to S3',
    schedule_interval='@daily',
)


def extract_data(**kwargs):
    session = SessionLocal()
    today = datetime.now().date()
    query = session.query(Molecule).filter(Molecule.created_at == today).all()  # Assuming 'created_at' column exists
    data = [{"identifier": mol.identifier, "smile": mol.smile} for mol in query]
    session.close()
    return data


def transform_data(ti, **kwargs):
    molecules = ti.xcom_pull(task_ids='extract_data')


    def calculate_properties(smile):
        mol = Chem.MolFromSmiles(smile)
        if not mol:
            return None
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        tpsa = rdMolDescriptors.CalcTPSA(mol)
        h_donors = rdMolDescriptors.CalcNumHBD(mol)
        h_acceptors = rdMolDescriptors.CalcNumHBA(mol)
        lipinski_pass = (mw <= 500 and logp <= 5 and h_donors <= 5 and h_acceptors <= 10)
        return [mw, logp, tpsa, h_donors, h_acceptors, lipinski_pass]

    transformed_data = []
    for mol in molecules:
        properties = calculate_properties(mol['smile'])
        if properties:
            transformed_data.append({
                "identifier": mol['identifier'],
                "smile": mol['smile'],
                "molecular_weight": properties[0],
                "logP": properties[1],
                "TPSA": properties[2],
                "H_donors": properties[3],
                "H_acceptors": properties[4],
                "lipinski_pass": properties[5]
            })

    df = pd.DataFrame(transformed_data)
    file_path = '/tmp/molecule_data.xlsx'
    df.to_excel(file_path, index=False)
    return file_path


def load_to_s3(ti, **kwargs):
    file_path = ti.xcom_pull(task_ids='transform_data')
    file_name = os.path.basename(file_path)
    s3_client.upload_file(file_path, S3_BUCKET, file_name)
    os.remove(file_path)

# Defined Task
extract_task = PythonOperator(
    task_id='extract_data',
    python_callable=extract_data,
    dag=dag,
)

transform_task = PythonOperator(
    task_id='transform_data',
    python_callable=transform_data,
    provide_context=True,
    dag=dag,
)

load_task = PythonOperator(
    task_id='load_to_s3',
    python_callable=load_to_s3,
    provide_context=True,
    dag=dag,
)

extract_task >> transform_task >> load_task
