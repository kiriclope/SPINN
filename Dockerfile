# Use an official GCC runtime as a parent image
FROM gcc:latest

# Set the working directory in the container
WORKDIR /usr/src/app

# Copy the current directory contents into the container
COPY . .

# Install any needed packages
RUN apt-get update && \
    apt-get install -y cmake make

# Compile the project
RUN cmake . && make

# Run the output program from the previous step
CMD ["./LifNet"]