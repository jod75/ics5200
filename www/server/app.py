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

class CustomEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, Decimal):
            return float(o)
        return super(CustomEncoder, self).default(o)

@main.route("/doExperiment/<molRegNo>", methods=["GET"])
def doExperiment(molRegNo):
    moleculeData = ics5200Engine.doExperiment(molRegNo)
    return json.dumps(moleculeData, cls=CustomEncoder)
    
@main.route("/smilesToSVG/<molRegNo>", methods=["GET"])
def smilesToSVG(molRegNo):    
    smiles = ics5200Engine.getSmiles(molRegNo)
    if len(smiles) > 0:
        logger.debug("smilesToSVG: " + smiles)
        return ChemInfo.smilesToSVG(smiles)
    else:
        return ""
        
@main.route("/testSmilesToSVG/<molRegNo>", methods=["GET"])
def testSmilesToSVG(molRegNo):    
    return ChemInfo.smilesToSVG(ics5200Engine.getTestSmiles(molRegNo))
        
@main.route("/getTestMolRegNos", methods=["GET"])
def getTestMolRegNos():
    return json.dumps(ics5200Engine.getTestMolRegs())
    
@main.route("/getTestTargets/<molRegNo>", methods=["GET"])
def getTestTargets(molRegNo):
    return json.dumps(ics5200Engine.getTestTargets(molRegNo))
    
@main.route("/getTestDataset/", methods=["GET"])
def getTestDataset():
    return json.dumps(ics5200Engine.getTestDataset())
    
@main.route("/doProteinExperiment/<compId>", methods=["GET"])
def doProteinExperiment(compId):
    proteinData = ics5200Engine.doProteinExperiment(compId)
    return json.dumps(proteinData, cls=CustomEncoder)

def createApp(sparkContext, datasetPath):
    global ics5200Engine

    ics5200Engine = ICS5200Engine(sparkContext, datasetPath)

    app = Flask(__name__)
    CORS(app)
    app.register_blueprint(main)
    return app

