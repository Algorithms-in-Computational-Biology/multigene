import StandardPCRPrimer as standard
import MultigenePCRPrimer as multigene

from wrapper import Multigene
from flask import Flask
from flask import render_template, request
from wtforms import Form, SelectField, TextAreaField, IntegerField, validators

app = Flask(__name__)

def get_sequences(inputs):
    sequences = []
    items = inputs.split('\r\n')
    for i in range(len(items)):
        if (i%2 == 1):
            sequences.append(items[i])
    return sequences


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
                 16, 
                 4, 
                 14, 
                 4, 
                 form.product_length.data]

        sequences = get_sequences(form.sequence.data)
        print(sequences)
        if (form.pcr_type.data == '1'):
            result = standard.design(sequences[0], param)
        
            return render_template('standard_result.html',result=result)
        else:
            #result = multigene.design(sequences, param)
            multigene = Multigene()
            result = multigene.design(sequences)

            return render_template('multigene_result.html',result=result)

    return render_template('parameter.html', form=form)


class ParameterForm(Form):
    pcr_type = SelectField('PCR type:', choices = [('1', 'Single sequence'), ('2', 'Multiple sequences')])
    sequence = TextAreaField('Paste your sequence (5\'&rarr;3\'):', [validators.length(max = 10000), 
                                                                     validators.Required("Sequence is required.")])

    min_length = IntegerField("min.", [validators.Required("Min. length is required.")], default=18)
    opt_length = IntegerField("opt.", [validators.Required("Opt. length is required.")], default=20)
    max_length = IntegerField("max.", [validators.Required("Max. length is required.")], default=24)

    min_content = IntegerField("min.", [validators.Required("Min. GC content is required.")], default=40)
    opt_content = IntegerField("opt.", [validators.Required("Opt. GC content is required.")], default=50)
    max_content = IntegerField("max.", [validators.Required("Max. GC content is required.")], default=60)

    min_tm = IntegerField("min.", [validators.Required("Min. melting temperature is required.")], default=50)
    opt_tm = IntegerField("opt.", [validators.Required("Opt. melting temperature is required.")], default=57)
    max_tm = IntegerField("max.", [validators.Required("Max. melting temperature is required.")], default=65)

    product_length = IntegerField('Product length:', [validators.Required("Product length is required.")], default=100)
