FROM sd2e/base:miniconda3-edge

WORKDIR /

ADD . /prot_stab_analysis_pipeline

# Start PEAR

RUN apt-get install -y build-essential autoconf automake libtool

RUN git clone https://github.com/stevschmid/PEAR

WORKDIR /PEAR

RUN ./autogen.sh
RUN ./configure
RUN make
RUN make install

RUN conda create --name pear_test

RUN bash -c "source activate pear_test && cd test && python test.py && source deactivate && rm -rf pear_test"

RUN which pear
# End PEAR

WORKDIR /

RUN conda env create -f /prot_stab_analysis_pipeline/environment.yml && rm -rf /opt/conda/pkgs/*

ENV PATH /opt/conda/envs/2018_prot_stab/bin:$PATH
ENV CONDA_DEFAULT_ENV 2018_prot_stab
ENV CONDA_PREFIX /opt/conda/envs/2018_prot_stab

CMD python /prot_stab_analysis_pipeline/scripts/compute_ec50_values_from_deep_sequencing_data.py