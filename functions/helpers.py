import os
import sys
import typing
import shutil
from pathlib import Path

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.text as mtext
import numpy as np
import math
import csv



def getProjectRoot() -> str:
    """Returns project root folder."""
    return os.path.normpath(str(Path(__file__).resolve().parents[1]))


def getRelativePath(target_path: str, base_path: str) -> str:
    """Get relative path from `base_path` to `target_path`. On Windows, if the `target_path` is not on the same drive as the `base_path`, return the `target_path`."""
    relpath = target_path
    try:
        relpath = os.path.relpath(target_path, base_path)
    except ValueError:
        pass
    return relpath


def applyBasicPlotStyle(plt: plt) -> plt:
    # Default font size to 18, font weight to medium, and font family to sans-serif Segoe UI
    plt.rcParams["font.size"] = 18
    plt.rcParams["font.weight"] = "medium"
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Segoe UI", "Myriad", "Tahoma", "DejaVu Sans", "Lucida Grande", "Verdana"]
    
    plt.rcParams["axes.titlesize"] = 22
    plt.rcParams["axes.titleweight"] = "bold"
    plt.rcParams["axes.labelsize"] = 19
    plt.rcParams["axes.labelweight"] = "medium"
    
    plt.rcParams["legend.fontsize"] = 14
    plt.rcParams["legend.title_fontsize"] = 16
    plt.rcParams["legend.framealpha"] = 1.0
    
    plt.rcParams["xtick.labelsize"] = 18
    plt.rcParams["ytick.labelsize"] = 18
    plt.rcParams["lines.linewidth"] = 2

    plt.rcParams["grid.alpha"] = 0.5
    plt.rcParams["grid.color"] = "black"

    return plt



class CurvedText(mtext.Text):
    """
    A text object that follows an arbitrary curve.
    """
    def __init__(self, x, y, text, axes, **kwargs):
        super(CurvedText, self).__init__(x[0],y[0],' ', **kwargs)

        axes.add_artist(self)

        ##saving the curve:
        self.__x = x
        self.__y = y
        self.__zorder = self.get_zorder()

        ##creating the text objects
        self.__Characters = []
        for c in text:
            if c == ' ':
                ##make this an invisible 'a':
                t = mtext.Text(0,0,'a')
                t.set_alpha(0.0)
            else:
                t = mtext.Text(0,0,c, **kwargs)

            #resetting unnecessary arguments
            t.set_ha('center')
            t.set_rotation(0)
            t.set_zorder(self.__zorder +1)

            self.__Characters.append((c,t))
            axes.add_artist(t)


    ##overloading some member functions, to assure correct functionality
    ##on update
    def set_zorder(self, zorder):
        super(CurvedText, self).set_zorder(zorder)
        self.__zorder = self.get_zorder()
        for c,t in self.__Characters:
            t.set_zorder(self.__zorder+1)

    def draw(self, renderer, *args, **kwargs):
        """
        Overload of the Text.draw() function. Do not do
        do any drawing, but update the positions and rotation
        angles of self.__Characters.
        """
        self.update_positions(renderer)

    def update_positions(self,renderer):
        """
        Update positions and rotations of the individual text elements.
        """

        #preparations

        ##determining the aspect ratio:
        ##from https://stackoverflow.com/a/42014041/2454357

        ##data limits
        xlim = self.axes.get_xlim()
        ylim = self.axes.get_ylim()
        ## Axis size on figure
        figW, figH = self.axes.get_figure().get_size_inches()
        ## Ratio of display units
        _, _, w, h = self.axes.get_position().bounds
        ##final aspect ratio
        aspect = ((figW * w)/(figH * h))*(ylim[1]-ylim[0])/(xlim[1]-xlim[0])

        #points of the curve in figure coordinates:
        x_fig,y_fig = (
            np.array(l) for l in zip(*self.axes.transData.transform([
            (i,j) for i,j in zip(self.__x,self.__y)
            ]))
        )

        #point distances in figure coordinates
        x_fig_dist = (x_fig[1:]-x_fig[:-1])
        y_fig_dist = (y_fig[1:]-y_fig[:-1])
        r_fig_dist = np.sqrt(x_fig_dist**2+y_fig_dist**2)

        #arc length in figure coordinates
        l_fig = np.insert(np.cumsum(r_fig_dist),0,0)

        #angles in figure coordinates
        rads = np.arctan2((y_fig[1:] - y_fig[:-1]),(x_fig[1:] - x_fig[:-1]))
        degs = np.rad2deg(rads)


        rel_pos = 10
        for c,t in self.__Characters:
            #finding the width of c:
            t.set_rotation(0)
            t.set_va('center')
            bbox1  = t.get_window_extent(renderer=renderer)
            w = bbox1.width
            h = bbox1.height

            #ignore all letters that don't fit:
            if rel_pos+w/2 > l_fig[-1]:
                t.set_alpha(0.0)
                rel_pos += w
                continue

            elif c != ' ':
                t.set_alpha(1.0)

            #finding the two data points between which the horizontal
            #center point of the character will be situated
            #left and right indices:
            il = np.where(rel_pos+w/2 >= l_fig)[0][-1]
            ir = np.where(rel_pos+w/2 <= l_fig)[0][0]

            #if we exactly hit a data point:
            if ir == il:
                ir += 1

            #how much of the letter width was needed to find il:
            used = l_fig[il]-rel_pos
            rel_pos = l_fig[il]

            #relative distance between il and ir where the center
            #of the character will be
            fraction = (w/2-used)/r_fig_dist[il]

            ##setting the character position in data coordinates:
            ##interpolate between the two points:
            x = self.__x[il]+fraction*(self.__x[ir]-self.__x[il])
            y = self.__y[il]+fraction*(self.__y[ir]-self.__y[il])

            #getting the offset when setting correct vertical alignment
            #in data coordinates
            t.set_va(self.get_va())
            bbox2  = t.get_window_extent(renderer=renderer)

            bbox1d = self.axes.transData.inverted().transform(bbox1)
            bbox2d = self.axes.transData.inverted().transform(bbox2)
            dr = np.array(bbox2d[0]-bbox1d[0])

            #the rotation/stretch matrix
            rad = rads[il]
            rot_mat = np.array([
                [math.cos(rad), math.sin(rad)*aspect],
                [-math.sin(rad)/aspect, math.cos(rad)]
            ])

            ##computing the offset vector of the rotated character
            drp = np.dot(dr,rot_mat)

            #setting final position and rotation:
            t.set_position(np.array([x,y])+drp)
            t.set_rotation(degs[il])

            t.set_va('center')
            t.set_ha('center')

            #updating rel_pos to right edge of character
            rel_pos += w-used


def compare_dictionaries(dict_1, dict_2, dict_1_name, dict_2_name, path=""):
    """Compare two dictionaries recursively to find non matching elements

    Parameters
    ----------
    - dict_1: **dict** 
        dictionary 1
    - dict_2: **dict**
        dictionary 2
    - dict_1_name: **str**
        name of dictionary 1
    - dict_2_name: **str**
        name of dictionary 2

    Returns
    ----------
    Information about the differences between the two dictionaries: **str**

    Example
    ----------
    
    ```python
    dict_1 = {'a': 1, 'b': 2, 'c': {'d': 3, 'e': 4}}
    dict_2 = {'a': 1, 'b': 2, 'c': {'d': 3, 'e': 5}}
    print(compare_dictionaries(dict_1, dict_2, 'dict_1', 'dict_2'))
    # Output: Value of dict_1[c][e] (4) not same as dict_2[c][e] (5)

    ```

    """
    err = ''
    key_err = ''
    value_err = ''
    old_path = path
    for k in dict_1.keys():
        path = old_path + "[%s]" % k
        if not dict_2.has_key(k):
            key_err += "Key %s%s not in %s\n" % (dict_1_name, path, dict_2_name)
        else:
            if isinstance(dict_1[k], dict) and isinstance(dict_2[k], dict):
                err += compare_dictionaries(dict_1[k],dict_2[k],'d1','d2', path)
            else:
                if dict_1[k] != dict_2[k]:
                    value_err += "Value of %s%s (%s) not same as %s%s (%s)\n"\
                        % (dict_1_name, path, dict_1[k], dict_2_name, path, dict_2[k])

    for k in dict_2.keys():
        path = old_path + "[%s]" % k
        if not dict_1.has_key(k):
            key_err += "Key %s%s not in %s\n" % (dict_2_name, path, dict_1_name)

    return key_err + value_err + err



def progress(value, length=40, title="", min_title_length=10, vmin=0.0, vmax=1.0, postfix="", auto_resize=True):
    """
    Text progress bar.
    
    Parameters
    ----------
    value : float
        Current value to be displayed as progress
    vmin : float
        Minimum value
    vmax : float
        Maximum value
    length: int
        Bar length (in character)
    title: string
        Text to be prepend to the bar
    postfix: string
        Text to be append at the end of the bar
    auto_resize: bool
        Auto adjust bar length to fit the screen size. If False, use the given length
    """

    LINE_UP = '\033[1A'
    LINE_CLEAR = '\x1b[2K'
    # Block progression is 1/8
    blocks = ["", "▏","▎","▍","▌","▋","▊","▉","█"]
    vmin = vmin or 0.0
    vmax = vmax or 1.0
    lsep, rsep = " ▏", "▕ "

    cols, _ = shutil.get_terminal_size(fallback = (length, 1))
    # limit cols to some % of terminal width
    cols = int(cols * 0.96)

    # remove title and/or postfix if the length is too long
    while len(title) + len(postfix) > cols - len(lsep) - len(rsep) - 10:
        if len(title) + len(postfix) <= min_title_length:
            break

        # remove postfix first, then if still too long, remove title gradually
        if len(postfix) > 0:
            postfix = ""
        else:
            title = title[:-1]

    # Auto adjust length to fit the screen by subtracting the length of title, lsep, rsep, percentage (10 char), and postfix
    if auto_resize:
        length = cols - len(title) - len(lsep) - len(rsep) - 10 - len(postfix)

    # Normalize value
    value = min(max(value, vmin), vmax)
    value = (value-vmin)/float(vmax-vmin)

    if length < 1:
        sys.stdout.write(LINE_UP + LINE_CLEAR + " (%.1f%%)" % (value*100) + "\n")
        sys.stdout.flush()
        return

    v = value*length
    x = math.floor(v) # integer part
    y = v - x         # fractional part
    base = 0.125      # 0.125 = 1/8
    prec = 3
    i = int(round(base*math.floor(float(y)/base),prec)/base)
    bar = "█"*x + blocks[i]
    n = length-len(bar)
    bar = lsep + bar + " "*n + rsep

    sys.stdout.write(LINE_UP + LINE_CLEAR + title + bar + postfix + " (%.1f%%)" % (value*100) + "\n")
    sys.stdout.flush()



def rotateCoordinates(coordinates: np.ndarray, rotation: typing.Union[int, float], flip_x: bool) -> np.ndarray:
    """
    Rotate the spatial coordinates of the cells.

    Parameters
    ----------
    - `coordinates` : **np.ndarray**
        Numpy array of the spatial coordinates of the cells.
    - `rotation` : **int** or **float**
        The degree of rotation to apply to the spatial coordinates.

    Returns
    -------
    - `coordinates` : **np.ndarray**
        Numpy array of the spatial coordinates of the cells after rotation.
    """
    rotation = math.radians(rotation)
    
    x = coordinates[:, 0]
    y = coordinates[:, 1]

    x_new = x * math.cos(rotation) - y * math.sin(rotation)
    y_new = x * math.sin(rotation) + y * math.cos(rotation)

    if flip_x:
        x_new = -x_new

    return np.column_stack((x_new, y_new))



def loadYaoMapMyCellsColors() -> pd.DataFrame:
    """
    Load ../ref_data/cluster_to_cluster_annotation_membership_color_20230830.csv and get the colors corresponding to all MapMyCells annotations.

    Returns
    -------
    - `colors_map` : **pd.DataFrame**
        A DataFrame containing the colors of all MapMyCells annotations.
    """
    colors_map = {}

    if not os.path.exists(os.path.join(getProjectRoot(), "ref_data", "cluster_to_cluster_annotation_membership_color_20230830.csv")):
        raise FileNotFoundError("File 'cluster_to_cluster_annotation_membership_color_20230830.csv' not found in 'ref_data' folder.")
    
    path = os.path.join(getProjectRoot(), "ref_data", "cluster_to_cluster_annotation_membership_color_20230830.csv")

    colors_map = pd.read_csv(path, index_col=0, header=0, sep=",")

    return colors_map


YAO_MAPMYCELLS_COLORS = loadYaoMapMyCellsColors()