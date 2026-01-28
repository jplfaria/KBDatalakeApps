FROM kbase/sdkpython:3.8.10
LABEL maintainer="chenry@anl.gov"

# -----------------------------------------
# Install system dependencies
# -----------------------------------------
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# -----------------------------------------
# Install KBUtilLib for shared utilities
# This provides common KBase functionality:
# - KBWSUtils: Workspace operations
# - KBGenomeUtils: Genome parsing and analysis
# - KBModelUtils: Metabolic model utilities
# - KBCallbackUtils: Callback server handling
# - SharedEnvUtils: Configuration and token management
# -----------------------------------------

# -----------------------------------------
# Copy module files
# -----------------------------------------
COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

RUN cd /kb/module && \
    git clone https://github.com/cshenry/KBUtilLib.git

WORKDIR /kb/module

# Compile the module
RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
