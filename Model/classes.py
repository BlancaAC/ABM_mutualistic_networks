#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 12:23:36 2024

@author: galeanojav
"""
from dataclasses import dataclass
import numpy as np

@dataclass
class Agent_Plant:
    id: int
    sp: str
    x: float
    y: float
    sp_complete: str

@dataclass
class Agent_Pol:
    id: int
    specie: str
    x: float
    y: float
    radioAccion: float

    def random_xy_pol(self, xmin, xmax, ymin, ymax):
        """" Move a pollinator using Brownian movements.
        
        Args:
            agent: The pollinator agent object.
            xmin: The minimum x-coordinate of the environment.
            xmax: The maximum x-coordinate of the environment.
            ymin: The minimum y-coordinate of the environment.
            ymax: The maximum y-coordinate of the environment.
        """
        m = self.radioAccion

        self.x = (self.x + np.random.uniform(-m, m)-xmin) % (xmax-xmin)+xmin
        self.y = (self.y + np.random.uniform(-m, m)-ymin) % (ymax-ymin)+ymin

class Environment_plant_v1:
    def __init__(self, df_plantPM, random_position=False, xmin=None, xmax=None, ymin=None, ymax=None):
        if random_position:
            # Update the x and y coordinates with random values within the provided limits.
            df_plantPM['X'] = np.random.uniform(xmin, xmax, len(df_plantPM))
            df_plantPM['Y'] = np.random.uniform(ymin, ymax, len(df_plantPM))
            
        self.plant_list = df_plantPM.apply(lambda row: Agent_Plant(row.Plant_id, row.Plant_sp, row.X, row.Y, row.Plant_sp_complete), axis=1).tolist()


class Environment_plant:
    def __init__(self, df_plantPM, random_position=False, regular_position=False, xmin=None, xmax=None, ymin=None, ymax=None):
        if random_position:
            # Update the x and y coordinates with random values within the provided limits.
            df_plantPM['X'] = np.random.uniform(xmin, xmax, len(df_plantPM))
            df_plantPM['Y'] = np.random.uniform(ymin, ymax, len(df_plantPM))
        elif regular_position:
            # Update the x and y coordinates with regular values within the provided limits.
            num_plantas = len(df_plantPM)
            num_filas = int(np.ceil(np.sqrt(num_plantas)))
            num_columnas = int(np.ceil(num_plantas / num_filas))
            x_vals = np.linspace(xmin, xmax, num_columnas)
            y_vals = np.linspace(ymin, ymax, num_filas)
            
            # Crear una lista de todas las combinaciones posibles de (x, y)
            xy_combinations = [(x, y) for y in y_vals for x in x_vals]
            
            # Asegurar que el número de combinaciones no exceda el número de plantas
            xy_combinations = xy_combinations[:num_plantas]
            
            for i in range(num_plantas):
                df_plantPM.iloc[i, df_plantPM.columns.get_loc('X')] = xy_combinations[i][0]
                df_plantPM.iloc[i, df_plantPM.columns.get_loc('Y')] = xy_combinations[i][1]
        
        self.plant_list = df_plantPM.apply(lambda row: Agent_Plant(row.Plant_id, row.Plant_sp, row.X, row.Y, row.Plant_sp_complete), axis=1).tolist()

class Environment_pol:
     def __init__(self, df_polPM):
         self.pol_list = df_polPM.apply(lambda row: Agent_Pol(row.Pol_id, row.Specie, row.x, row.y, row.Radius), axis=1).tolist()

