import boto3
import json

# Initialize a session using your credentials
session = boto3.Session(
    aws_access_key_id='<YOUR_ACCESS_KEY>',
    aws_secret_access_key='<YOUR_SECRET_KEY>',
    region_name='<YOUR_REGION>'
)

# Initialize the Lambda client
lambda_client = session.client('lambda')

# Define the payload (input event)
payload = {
    "names": ["Alice", "Bob", "Charlie"]
}


# Invoke the Lambda function
response = lambda_client.invoke(
    FunctionName='HelloStudentFunction',
    InvocationType='RequestResponse',
    Payload=json.dumps(payload)
)

# Read the response
response_payload = json.load(response['Payload'])
print(response_payload)
