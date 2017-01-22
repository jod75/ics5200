#!/usr/bin/env python

##################################################################################
# app.py
from flask import Blueprint
main = Blueprint('main', __name__)

import json
from engine import ICS5200Engine

import logging
logging.basicConfig(filename = "~/ics5200.log", level=logging.INFO)
logger = logging.getLogger(__name__)

from flask import Flask, request
from flask_cors import CORS, cross_origin

from decimal import Decimal

from cheminfo import ChemInfo
from moleculehelper import *

class CustomEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, Decimal):
            return float(o)
        return super(CustomEncoder, self).default(o)

@main.route("/getTestLigandsDS/", methods=["GET"])
def getTestLigandsDS():
    return json.dumps(ics5200Engine.getTestLigandsDS(), cls=CustomEncoder)

@main.route("/getTestProteinsDS/", methods=["GET"])
def getTestProteinsDS():
    return json.dumps(ics5200Engine.getTestProteinsDS(), cls=CustomEncoder)

@main.route("/getSmilesSVG/<molRegNo>", methods=["GET"])
def getSmilesSVG(molRegNo):    
    smiles = ics5200Engine.getSmiles(molRegNo)
    if len(smiles) > 0:
        logger.debug("smilesToSVG: " + smiles)
        return ChemInfo.smilesToSVG(smiles)
    else:
        return ""

@main.route("/getLigandTestTargets/<molRegNo>", methods=["GET"])
def getLigandTestTargets(molRegNo):
    return json.dumps(ics5200Engine.getLigandsTestTargets(molRegNo))

@main.route("/doLigandExperiment/<molRegNo>/<fingerprint>/<similarity>/<threshold>", methods=["GET"])
def doLigandExperiment(molRegNo, fingerprint, similarity, threshold):
    return json.dumps(ics5200Engine.doLigandExperiment(molRegNo, LigandHelper, fingerprint, similarity, float(threshold)), cls=CustomEncoder)

@main.route("/getProteinTestBindings/<compId>", methods=["GET"])
def getProteinTestBindings(compId):
    return json.dumps(ics5200Engine.getProteinTestBindings(compId), cls=CustomEncoder)

@main.route("/doProteinExperiment/<compId>", methods=["GET"])
def doProteinExperiment(compId):
    return json.dumps(ics5200Engine.doProteinExperiment(compId), cls=CustomEncoder)

def createApp(sparkContext, datasetPath):
    global ics5200Engine

    ics5200Engine = ICS5200Engine(sparkContext, datasetPath)

    app = Flask(__name__)
    CORS(app)
    app.register_blueprint(main)
    return app

