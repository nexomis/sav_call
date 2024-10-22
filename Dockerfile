FROM ubuntu:22.04 AS build

# Install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    g++ \
    libhts-dev \
    make \
    && rm -rf /var/lib/apt/lists/* \
    mkdir /app

# Clone the repository
WORKDIR /app
COPY src/ /app/src
COPY Makefile /app/Makefile
RUN  make

# Stage 2: Final
FROM ubuntu:22.04

# Install runtime dependencies
RUN export DEBIAN_FRONTEND=noninteractive \ 
  && apt-get update \
  && apt-get -y install --no-install-recommends \
    libstdc++6 \
    libc6 \
    libhts3 \
    && rm -rf /var/lib/apt/lists/*

# Copy the built executable from the build stage
COPY --from=build /app/sav_call /usr/local/bin/sav_call
