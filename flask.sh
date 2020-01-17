#!/bin/bash

cd ./gui/Flask/
export FLASK_APP=app.py
export FLASK_ENV=development
. venv/bin/activate
python -m flask run
