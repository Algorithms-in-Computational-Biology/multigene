from wrapper import Multigene
from flask import Flask
from flask import render_template, request
from wtforms import Form, TextAreaField, IntegerField, DecimalField, validators

app = Flask(__name__)

def get_sequences(inputs):
    sequences = []
    items = inputs.split('\r\n')
    for i in range(len(items)):
        if (i%2 == 1):
            sequences.append(items[i])
    return sequences

def convert(items):
    l = []
    item = items
    while(item):
        l.append(item)
        item = item.contents.next
    return l

@app.route('/primer-design', methods=['GET', 'POST'])
def search():
    form = ParameterForm(request.form) 
    if (request.method == 'POST' and form.validate()):
        param = [form.min_length.data, 
                 form.max_length.data, 
                 form.min_content.data, 
                 form.max_content.data, 
                 form.min_tm.data, 
                 form.max_tm.data, 
                 form.min_product_length.data,
                 form.max_product_length.data,
                 form.primer_concentration.data,
                 form.template_concentration.data]

        multigene = Multigene()
        sequences = get_sequences(form.sequence.data)
        #result = multigene.design(sequences, param)
        result = convert(multigene.design(sequences))

        return render_template('multigene_result.html', result=result)

    return render_template('parameter.html', form=form)


class ParameterForm(Form):
    sequence = TextAreaField('Sequences (5\'&rarr;3\'): ', [validators.length(max = 10000), validators.Required("Sequences are required.")])

    min_length = IntegerField("Min.", [validators.Required("Minimum length is required.")], default=18)
    max_length = IntegerField("Max.", [validators.Required("Maximum length is required.")], default=24)

    min_content = IntegerField("Min.", [validators.Required("Minimum GC content is required.")], default=40)
    max_content = IntegerField("Max.", [validators.Required("Maximum GC content is required.")], default=60)

    min_tm = IntegerField("Min.", [validators.Required("Minimum melting temperature is required.")], default=50)
    max_tm = IntegerField("Max.", [validators.Required("Maximum melting temperature is required.")], default=65)

    min_product_length = IntegerField('Min.', [validators.Required("Minimum product length is required.")], default=100)
    max_product_length = IntegerField('Max.', [validators.Required("Maximum product length is required.")], default=1200)

    primer_concentration = DecimalField('Primer', [validators.Required("Primer concentration is required.")], places=6, default=0.000001)
    template_concentration = DecimalField('Template', [validators.Required("Template concentration is required.")], places=6, default=0.000001)
