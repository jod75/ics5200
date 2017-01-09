#!/usr/bin/env python

##################################################################################
# app.py
from flask import Blueprint
main = Blueprint('main', __name__)

import json
from engine import ICS5200Engine

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

from flask import Flask, request

@main.route("/<int:userId>/sample/<int:moleculeId>", methods=["GET"])
def sample(userId, moleculeId):
    logger.debug("User %s requested data for molecule%s", userId, moleculeId)
    moleculeData =ics5200Engine.getMoleculeInfo(userId, [moleculeId])
    return json.dumps(moleculeData)

def createApp(sparkContext, datasetPath):
    global ics5200Engine

    ics5200Engine = ICS5200Engine(sparkContext, datasetPath)

    app = Flask(__name__)
    app.register_blueprint(main)
    return app

