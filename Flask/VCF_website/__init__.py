
from flask import Flask, Blueprint
from flask_sqlalchemy import SQLAlchemy
from flask_session import Session



app = Flask(__name__)

app.config['SECRET_KEY'] = '718414b2da90be33fcfdf92d803e5dc6718414b2da90be33fcfdf92d803e5dc6'
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///site.db'
db = SQLAlchemy(app)



app.config['SESSION_TYPE'] ='filesystem'
sess = Session()
sess.init_app(app)

from VCF_website.main.routes import main
from VCF_website.query_results.routes import query_results
from VCF_website.statistics.routes import statistics
app.register_blueprint(main)
app.register_blueprint(query_results)
app.register_blueprint(statistics)


