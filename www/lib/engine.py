#!/usr/bin/env python

##################################################################################
# engine.py
import os
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ICS5200Engine:
    def getMoleculeInfo(self, userId, moleculeId):
        #requested_movies_RDD = self.sc.parallelize(movie_ids).map(lambda x: (user_id, x))
        # Get predicted ratings
        #ratings = self.__predict_ratings(requested_movies_RDD).collect()

        return "Hello"

    def __init__(self, sc, dataset_path):
        """Init the recommendation engine given a Spark context and a dataset path
        """

        logger.info("Starting up the ICS5200 Engine: ")

        self.sc = sc

        # Load ratings data for later use
        logger.info("Loading data...")

