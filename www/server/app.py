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

@main.route("/<int:userId>/sample/<int:moleculeId>", methods=["GET"])
def sample(userId, moleculeId):
    logger.debug("User %s requested data for molecule%s", userId, moleculeId)
    moleculeData =ics5200Engine.getMoleculeInfo(userId, moleculeId)
    return json.dumps(moleculeData, cls=CustomEncoder)
    
@main.route("/smilesToSVG/<molregno>", methods=["GET"])
def smilesToSVG(molregno):    
    smiles = ics5200Engine.getSmiles(molregno)
    if len(smiles) > 0:
        logger.debug("smilesToSVG: " + smiles[0][1])
        return ChemInfo.smilesToSVG(smiles[0][1])
    else:
        return ""
        
@main.route("/testSmilesToSVG/<int:molIndex>", methods=["GET"])
def testSmilesToSVG(molIndex):    
    return ChemInfo.smilesToSVG(ics5200Engine.getTestSmiles(molIndex))
        
@main.route("/getTestMolRegNos", methods=["GET"])
def getTestMolRegNos():
    return json.dumps(ics5200Engine.getTestMolRegNos())

def createApp(sparkContext, datasetPath):
    global ics5200Engine

    ics5200Engine = ICS5200Engine(sparkContext, datasetPath)

    app = Flask(__name__)
    CORS(app)
    app.register_blueprint(main)
    return app

