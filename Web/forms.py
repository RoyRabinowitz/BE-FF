# AUTHOR: Shiri Almog , shirialmog1@gmail.com
from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileAllowed
from wtforms import StringField, SubmitField, SelectField, BooleanField,FormField
from wtforms.validators import DataRequired, Length

class rsIDform(FlaskForm):
    rsID = StringField("If you know the rsID, enter it here:",id='rsID',validators=[DataRequired()])
    submitrsID = SubmitField('Submit')

class DNAform(FlaskForm):
    rsID=StringField("If you know the rsID, enter it here:")
    PAM=StringField("PAM",id='PAM',validators=None)
    PAMstart = StringField("Activity window (distance from PAM)",id='start')
    PAMend = StringField("end",id='end')
    pambase=SelectField('BE type',choices=[('0','Select'),('1','CBE (C to T)'),('2','ABE (A to G)')],default='0',id="fromto")
    pamstream=SelectField('PAM orientation',choices=[('0','Select'),('U','Upstream'),('D','Downstream')],default='0',id="stream")

    upSeq = StringField("Upstream sequence:",id="Useq" ,validators=[DataRequired(), Length(min=25)])
    downSeq = StringField("Downstream sequence:", id="Dseq",
                          validators=[DataRequired(), Length(min=25)])
    mutation = StringField("Variant:",id="Mut", validators=[DataRequired(), Length(min=1, max=1)])
    WT = StringField("Reference:", id="Var",validators=[DataRequired(), Length(min=1, max=1)])
    readingFrame = StringField('Reading Frame:', id="RF", validators=None)
    #readingFrame = StringField('Reading Frame:',id="RF",validators=[DataRequired(),Length(min=1,max=1)])
    #showReadingFrame = BooleanField("Reading Frame:", default=True,validators=None,id='isRF')
    submit=SubmitField('Enter')
    #submitrsID = SubmitField('Submit')

class beffForm(FlaskForm):
    rsIDf=FormField(rsIDform)
    DNAf=FormField(DNAform)

class contactForm(FlaskForm):
    message = StringField("Leave us a message here:", validators=[DataRequired()])
    submit = SubmitField('Enter')

class crisPAMform(FlaskForm):
    referencePAM = StringField("Enter reference Sequence here:", validators=[DataRequired(),Length(min=25)])
    variantPAM = StringField("Enter variant Sequence here:", validators=[DataRequired(),Length(min=25)])
    submit = SubmitField('Enter')

class uploadForm(FlaskForm):
    inputFile=FileField('Upload file here', validators=[DataRequired(), FileAllowed(['csv'])])
    submit=SubmitField('Upload')