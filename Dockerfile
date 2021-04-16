# This Dockerfile uses multi-stage building to minimize size
# while maintaining all necessary dependencies. The build
# stage brings in dependencies necessary for building setsm
# and then compiles the source code. The run stage brings the
# executable created from the build stage into a new image and
# includes required dependencies

# Command Line arguments
ARG VERSION=16.04
ARG COMPILER=gnu

# ------------- BUILD STAGE ------------- #
FROM ubuntu:$VERSION as builder

# Bring arguments in from command line
ARG COMPILER
ARG VERSION

# Specify that the bash shell will be used
SHELL ["/bin/bash", "-c"]

# Bring in dependencies
RUN apt-get update && apt-get install --no-install-recommends -y \
    libgeotiff-dev \
    libgeotiff[0-9]+ \
    g++ \
    git \
    ca-certificates \
    make \
    wget \
    gnupg2 \
    apt-transport-https \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install libproj-dev
RUN apt-get update && apt-get install -y libproj-dev

# Create file that holds compiler specific paths
RUN touch /opt/compilerpath

# If building Intel version, then install Intel compiler
RUN if [ "$COMPILER" = 'intel' ]; then \
wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB; \
apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB; \
rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB; \
echo "deb https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list; \
apt-get update; \
apt-get install -y intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic; \
source /opt/intel/oneapi/setvars.sh; \
echo $PATH > /opt/compilerpath; \
apt-cache pkgnames intel | grep libiomp; \
echo **DONE**; \
fi

# Change working directory to /opt
WORKDIR /opt

# Copy all files into /opt
COPY ./* /opt/

# Update path and compiler
ENV PATH="/opt:${PATH}"
ENV COMPILER=$COMPILER

# Bring in compiler-specific paths, and then make the setsm executable
RUN PATH="$(cat compilerpath):$PATH"; make COMPILER=$COMPILER INCS=-I/usr/include/geotiff

# ------------- RUN STAGE ------------- #
FROM ubuntu:$VERSION as runner

# Bring in dependencies
RUN apt-get update && apt-get install --no-install-recommends -y \
        libgeotiff[0-9]+ \
        libgomp1 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Make setsmdir and bring files into it
RUN mkdir /opt/setsmdir
COPY --from=builder /opt/setsm /opt/setsmdir/

# Update path
ENV PATH="/opt/setsmdir:${PATH}"

# When image is run, setsm is executed
CMD setsm
