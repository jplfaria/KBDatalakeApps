FROM python:3.10-slim-bullseye
LABEL us.kbase.python="3.10"
LABEL maintainer="chenry@anl.gov"

ENV TZ=Etc/UTC
ENV DEBIAN_FRONTEND=noninteractive
ENV PIP_PROGRESS_BAR=off

# -----------------------------------------
# Install system dependencies
# -----------------------------------------
RUN apt-get update
RUN apt-get install -y build-essential \
                       git \
                       openjdk-11-jre \
                       unzip \
                       htop \
                       wget \
                       curl \
                       gcc \
                       cmake

#RUN rm -rf /var/lib/apt/lists/*

# Install uv (goes to /root/.local/bin by default)
RUN wget -qO- https://astral.sh/uv/install.sh | sh
RUN mkdir -p /opt/env
RUN /root/.local/bin/uv venv --python 3.10 /opt/env/berdl_genomes

# setup modelseedpy ml env
RUN /root/.local/bin/uv venv --python 3.10 /opt/env/modelseedpy_ml
RUN /root/.local/bin/uv pip install --python /opt/env/modelseedpy_ml --no-progress \
    git+https://github.com/ModelSEED/ModelSEEDpy.git@d109a897ac4b22c8ae2a78c963ad6ce01af9ce6a
RUN mkdir -p /kb/module/data
RUN wget -O /kb/module/data/knn_ACNP_RAST_full_01_17_2023_features.json -q \
    https://bioseed.mcs.anl.gov/~fliu/modelseedpy/knn_ACNP_RAST_full_01_17_2023_features.json
RUN wget -O /kb/module/data/knn_ACNP_RAST_full_01_17_2023.pickle -q \
    https://bioseed.mcs.anl.gov/~fliu/modelseedpy/knn_ACNP_RAST_full_01_17_2023.pickle

# Copy in the SDK
COPY --from=kbase/kb-sdk:1.2.1 /src /sdk
RUN sed -i 's|/src|/sdk|g' /sdk/bin/*
ENV PATH=/sdk/bin:$PATH

ADD requirements_kbase.txt /tmp/requirements_kbase.txt
RUN /usr/local/bin/pip install -r /tmp/requirements_kbase.txt
ADD biokbase /opt/conda/lib/python3.11/site-packages
ADD biokbase/user-env.sh /kb/deployment/user-env.sh
# add rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

RUN git clone https://github.com/soedinglab/MMseqs2.git /opt/MMseqs2
RUN git clone https://github.com/bluenote-1577/skani.git /opt/skani

ENV RUSTUP_INIT_SKIP_PATH_CHECK=yes
ENV CARGO_NET_GIT_FETCH_WITH_CLI=true
ENV PATH="/root/.cargo/bin:${PATH}"
RUN rustc --version && cargo --version
WORKDIR /opt/skani
RUN cargo install --path . --root ~/.cargo

RUN mkdir -p /deps
RUN echo '0' >/dev/null && cd /deps && \
	git clone https://github.com/ModelSEED/ModelSEEDDatabase.git && \
    cd ModelSEEDDatabase && git checkout 3346b71a34bc9d8c5a365b71d5a2959ffbe6c26e
RUN pip install --upgrade pip
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

RUN wget -O /tmp/mmseqs_avx2.tar.gz -q https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz
RUN mkdir -p /opt/mmseqs_avx
RUN tar -xvf /tmp/mmseqs_avx2.tar.gz -C /opt/mmseqs_avx
ENV PATH=/opt/mmseqs_avx/mmseqs/bin:$PATH

ADD requirements.txt /tmp/requirements.txt

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module
RUN /usr/local/bin/pip install -r /tmp/requirements.txt

RUN /root/.local/bin/uv pip install --python /opt/env/berdl_genomes --no-progress -r /kb/module/berdl/requirements.txt

# @chenry
WORKDIR /deps
RUN echo 'run this'
RUN git clone https://github.com/cshenry/ModelSEEDpy.git && pip install --use-deprecated=legacy-resolver -e ModelSEEDpy
RUN echo '0' >/dev/null && cd /deps && \
    git clone https://github.com/cshenry/cobrakbase.git && \
    cd cobrakbase && git checkout 68444e46fe3b68482da80798642461af2605e349
RUN echo '6' >/dev/null && cd /deps && \
    git clone https://github.com/cshenry/KBUtilLib.git
RUN echo '0' >/dev/null && cd /deps && \
    git clone https://github.com/kbaseapps/cb_annotation_ontology_api.git

WORKDIR /kb/module

# Compile the module
RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
