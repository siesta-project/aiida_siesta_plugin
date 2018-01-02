FROM aiida_siesta_code:0.10.1

COPY . /code/aiida_siesta_plugin
WORKDIR /code/aiida_siesta_plugin

RUN pip install -r ./test_requirements.txt \
    && pip install -e .
