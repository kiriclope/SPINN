# Start from a base Image
FROM ubuntu:latest

# Install essential libraries
RUN apt-get update && apt-get install -y \
    g++ \
    cmake \
    make

# Set the working directory to /app
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY . /app

# Build the app
RUN cmake .
RUN cmake --build .

# Run the output app using the command line
CMD ["./LifNet"]
