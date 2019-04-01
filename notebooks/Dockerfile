FROM jupyter/scipy-notebook:50247fea877a

WORKDIR /mnt/code

RUN conda install -y -c rdkit rdkit

COPY requirements.txt requirements.txt
RUN pip install --upgrade pip
RUN pip install -r requirements.txt

RUN ipython3 profile create
RUN echo 'import pandas as pd; pd.set_option("display.width", 160)' >> /home/jovyan/.ipython/profile_default/ipython_config.py

RUN jupyter contrib nbextension install --user
RUN jupyter nbextensions_configurator enable --user
RUN jupyter nbextension enable --py widgetsnbextension
RUN jupyter nbextension enable execute_time/ExecuteTime
RUN jupyter nbextension enable codefolding/main
RUN jupyter nbextension enable codefolding/edit
RUN jupyter nbextension enable toc2/main
RUN jupyter nbextension enable collapsible_headings/main
