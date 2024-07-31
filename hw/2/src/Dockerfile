# Use Miniconda base image
FROM continuumio/miniconda3

# Set working directory
WORKDIR /app

# Install RDKit
RUN conda install conda-forge::rdkit

# Install dependencies
COPY requirements.txt .
RUN pip3 install -r requirements.txt

# Activate environment and set it as default
RUN echo "conda activate myenv" >> ~/.bashrc
ENV PATH /opt/conda/envs/myenv/bin:$PATH

# Copy application code
COPY . .

# Expose the port FastAPI will run on
EXPOSE 8000

# Run the FastAPI application
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
