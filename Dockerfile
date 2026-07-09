# Use a stable, official Python core image
FROM python:3.12-slim

# Install core system packages required by Snakemake and SciPy compiling
RUN apt-get update && apt-get install -y \
    build-essential \
    graphviz \
    git \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory inside the container
WORKDIR /workspace

# Copy requirements over first (helps speed up builds via caching)
COPY requirements.txt .
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# Copy the rest of your project directories (lcr_webapp, scripts, etc.)
COPY . .

# Expose the default network port Render uses
EXPOSE 10000

# Start your Flask app via Gunicorn bound to the container port
CMD ["gunicorn", "--bind", "0.0.0.0:10000", "lcr_webapp.app:app"]