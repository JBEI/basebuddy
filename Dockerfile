FROM python:3.8
EXPOSE 8080
WORKDIR /app
COPY requirements.txt .
ADD streamlit-searchbox /app/streamlit-searchbox
RUN pip install -r requirements.txt
COPY . ./
RUN cd data && tar -xvf cocoput_table.tsv.tar.gz && cd ..
ENTRYPOINT ["streamlit", "run", "streamlit.py", "--server.port=8080", "--server.address=0.0.0.0"]
