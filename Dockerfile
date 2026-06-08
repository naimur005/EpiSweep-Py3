FROM python:3.10-bullseye

RUN apt-get update -qq 
RUN apt-get install -y -qq --no-install-recommends r-base && apt-get clean
RUN apt-get update && apt-get install -y ncbi-blast+
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

# Set environment variables
ENV PYTHONPATH="/app"
ENV EPISWEEP="/app"

# Copy the CPLEX installer and response file to the Docker image
COPY cplex_studio2211.linux_x86_64.bin /tmp/cplex_installer.bin
COPY cplex_response_file.properties /tmp/cplex_response_file.properties

# Make the installer executable and run it in silent mode using the response file
RUN chmod +x /tmp/cplex_installer.bin && \
    /tmp/cplex_installer.bin -f /tmp/cplex_response_file.properties && \
    rm /tmp/cplex_installer.bin /tmp/cplex_response_file.properties

# Set the Path & install the CPLEX Python API for Python 3.10
ENV PATH="/opt/ibm/ILOG/CPLEX_Studio2211/cplex/bin/x86-64_linux/:$PATH"
RUN cd /opt/ibm/ILOG/CPLEX_Studio2211/cplex/python/3.10/x86-64_linux/ && \
    python3 setup.py install

# Download and install JDK 17
COPY jdk-17_linux-x64_bin.tar.gz /tmp/
RUN mkdir -p /usr/lib/jvm && \
    tar -xzf /tmp/jdk-17_linux-x64_bin.tar.gz -C /usr/lib/jvm && \
    rm /tmp/jdk-17_linux-x64_bin.tar.gz
ENV JAVA_HOME=/usr/lib/jvm/jdk-17.0.10
ENV PATH="$JAVA_HOME/bin:$PATH"

# Copy OSPREY files and set permissions
COPY misc/OSPREY3.2 /app/misc/OSPREY3.2
ENV OSPREY_PATH="/app/misc/OSPREY3.2"
RUN chmod +x $OSPREY_PATH/*

RUN R --version
RUN python --version
RUN java -version

# Set the working directory in the container
WORKDIR /app

# Copy the current directory contents into the container
COPY . .

# Install any Python dependencies
RUN pip install --no-cache-dir -r requirements.txt
RUN pip install --upgrade pip

# Keep the container running 
CMD ["tail", "-f", "/dev/null"]
