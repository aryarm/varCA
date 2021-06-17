FROM snakemake/snakemake:v5.24.2

# work in the home directory
WORKDIR /home

# copy all the files to the container
COPY . .

# install the pipeline's dependencies
RUN snakemake --use-conda --conda-create-envs-only -j

# run varCA
# specify as a CMD, so it can be overridden by the user
CMD ["./run.bash"]
