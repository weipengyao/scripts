__all__ = ["dataAnalysisUtil"]

#A class that comes with functions that make data analysis easier
#It plots, returns arrays and what have you

###############################################################

import os, sys, string, re, time, shutil, types, glob, socket, math, getopt
import yt
import numpy as np
import matplotlib.pyplot as plt


import os.path



#Figure out why docstrings are being weird on other machines
#consistent parameters
#more detailed documentation
PLog = True
#PLog just tells all functions whether they should print what they're doing. It's on by default, but you can pass false to init to have fewer prints.
#You can also specify it in the function itself as an argument if you want to overwrite the one provided in init.
def init(PrintStat = True):
    global PLog
    PLog = PrintStat

    
#Loads up a given file using YT and returns it    
def loadup(fileName):
    """
    Loads and returns a given file in YT.
    You must load up an file before you use ANY function that requires a "ds" object. 
    """
    ds = yt.frontends.flash.data_structures.FLASHDataset(filename=fileName)
    #ds = yt.frontends.flash.data_structures.FLASHDataset(filename=fileName)
    print("Loaded up %s" % fileName)
    return ds

#turns a given object into a numpy array
def toNumpyArray(frb,field,location='flash',POutput=None):
    """
    A function that quickly converts an FRB to a numpy array.
    
    Parameters:
        frb: The fixed resolution buffer to convert, although it can theoretically be used to convert non-FRB objects of a similar form to numpy arrays.
        field: The field to pull out of the FRB. It defulats to pulling from flash, but can be specified as an argument.
        location: Where the field is located, by default it goes to flash, but you can change it to "io" or "gas" for example.
    Returns:
        A numpy array.
    """
    npArray = np.array(frb[(location,field)])
    if ((POutput == None) and (PLog == True)) or (POutput == True): print("Generating numpy array")
    return npArray

#outputs the geometry of the given ds. 2D cylinder, 3D Cartesian. mostly for gathering information at runtime
def fetchGeometry(ds,POutput=None):
    """
    A function that returns the geometry of a dataset.
    
    Parameters:
        ds: The dataset object to find the geometry of.

    Returns:
        The geometry type as a string.
    """
    geoType = ds.geometry
    geoDimensionality = ds.dimensionality
    geometry = str(geoDimensionality) + "D " + str(geoType)
    geometry = geometry.upper()
    if ((POutput == None) and (PLog == True)) or (POutput == True): print("Geometry found to be %s" % geometry)
    return geometry
    
#finds the number of cells at a given refinement level. If includeZ is true it will include the z axis in the list too
#although this isn't super useful, it's more information that the user may want
#remove include z
def findNumCells(ds, ref=0,POutput=None):
    """
    A function that finds the number of cells at a given refinement level.
    
    Parameters:
        ds: The dataset object to count.
        ref: Refinement level to use when counting, if the provided refinment level is to high it will default to the maximum possible level.
        POutput: Whether or not it should print runtime info, on by default as defined by init().
    Returns:
        A list containing the needed sizes. comes in the form [X] for 1D [X,Y] for 2D and [X,Y,Z] for 3D.
    """
    tref = checkPossibleLevel(ds,ref)
    
    
    nb_pts_x = ds.domain_dimensions[0]*2**(tref)
    
    if ds.dimensionality == 1:
       if ((POutput == None) and (PLog == True)) or (POutput == True):
        print("Number of cells (X): ", [nb_pts_x])
       return [nb_pts_x]
    elif ds.dimensionality == 2:
       nb_pts_y = ds.domain_dimensions[1]*2**(tref)
       if ((POutput == None) and (PLog == True)) or (POutput == True):
        print("Number of cells (X,Y): ", [nb_pts_x,nb_pts_y])
       return [nb_pts_x,nb_pts_y] 
    else:
       nb_pts_y = ds.domain_dimensions[1]*2**(tref)
       nb_pts_z = ds.domain_dimensions[2]*2**(tref)
       if ((POutput == None) and (PLog == True)) or (POutput == True):
        print("Number of cells (X,Y,Z): ", [nb_pts_x,nb_pts_y,nb_pts_z]) 
       return [nb_pts_x,nb_pts_y,nb_pts_z] 
        
    
#returns a list with the edges of the container on x and y by default, although it will include z as well if told to so
#it should just return all of them by defauly (xx for 1d, xx,yy for 2d, ect)
def findFullSize(ds,POutput=None):    
    """
    A function that returns a list containing the edges of the dataset object
    
    Parameters:
        ds: The dataset object to operate on, the return changes size and contents depending on the dimensionality of the data.
        POutput: Whether or not it should print runtime info, on by default as defined by init().
    Returns:
        A list containing the needed sizes, will return (xleft,xright) for 1D, (xleft,xright,yleft,yright) for 2D and (xleft,xright,yleft,yright,zleft,zright) for 3D 
    """
    if ds.dimensionality == 1:
        if ((POutput == None) and (PLog == True)) or (POutput == True):
            print("Domain left and right edge(X): ", (ds.domain_left_edge[0], ds.domain_right_edge[0]))
        return (ds.domain_left_edge[0], ds.domain_right_edge[0])
    elif ds.dimensionality == 2:
        if ((POutput == None) and (PLog == True)) or (POutput == True):
            print("Domain left and right edge(X): ", (ds.domain_left_edge[0], ds.domain_right_edge[0]))
            print("Domain left and right edge(Y): ", (ds.domain_left_edge[1], ds.domain_right_edge[1]))
        return (ds.domain_left_edge[0], ds.domain_right_edge[0], ds.domain_left_edge[1], ds.domain_right_edge[1])
    else:
        if ((POutput == None) and (PLog == True)) or (POutput == True):
            print("Domain left and right edge(X): ", (ds.domain_left_edge[0], ds.domain_right_edge[0]))
            print("Domain left and right edge(Y): ", (ds.domain_left_edge[1], ds.domain_right_edge[1]))
            print("Domain left and right edge(Z): ", (ds.domain_left_edge[2], ds.domain_right_edge[2]))
        return (ds.domain_left_edge[0], ds.domain_right_edge[0], ds.domain_left_edge[1], ds.domain_right_edge[1], ds.domain_left_edge[2],ds.domain_right_edge[2])

      
#Slices data along a given axis at given coordinates      
def sliceData(ds, axis, coords, center=None, fieldParameters=None, source=None,POutput=None):
    """
    A function that slices a given dataset object.
    NOTE: 2D cylinders can ONLY be cut along axis 2, anything else will raise an error
    
    Parameters:
        ds: The dataset object to slice.
        axis: The axis to slice along. If the given axis is not possible, it will raise an error. Can be 0,1,2 for x,y,z. On Cylindrical data you can only like along axis 2
        coords: The coordinate where the slice will be made.
        center: An optional parameter that passes the center to fields that use it.
        fieldParameters: A dictionary of field parameters than can be accessed by derived fields.
        source: Draw the selection from the provided data source rather than all data associated with the dataset.
        POutput: Whether or not it should print runtime info, on by default as defined by init().
    Returns:
        A slice object
    """
    if axis > 2:
        raise ValueError("impossible slice")
    if (ds.dimensionality == 2) and (axis != 2):
        raise TypeError("Cannot slice 2D data along any axis except 2")
    if ds.dimensionality == 1:
        raise TypeError("Cannot slice 1D data")
    if (ds.domain_right_edge[axis] < coords) or (ds.domain_left_edge[axis] > coords):
        raise TypeError("Coordinate out of bounds")
    if (center != None):
        if checkPossibleLocation(ds,center,POutput) == False:
            raise TypeError("Coordinates for center are out of bounds")                   
    slc = ds.slice(axis, coords, center=center, field_parameters=fieldParameters, data_source=source)
    if ((POutput == None) and (PLog == True)) or (POutput == True): print("Generating slice along axis %s " % axis)
    return slc
    
    
    #test FRB resLevel order?
def makeFRB(source, bounds, level=0, antialias=True, periodic=False,POutput=None):
    """
    A function that turns a given 2D object into a fixed resolution buffer
    
    Parameters:
        source: The 2D object to turn into an FRB, it can take a slice a projection or a covering grid. CHECK FOR COVERING GRID
        bounds: Bounds are the min and max in the image plane that we want our image to cover.  It's in the order of (xmin, xmax, ymin, ymax)?, where the coordinates are all in the appropriate code units. It can also automatically use the full domain with the strings "XY", "YZ" "XZ"
        resLevel: The size of the image to generate. It will take a given refinement level and use that to find the buffersize with that 
        antialias: Determines whether or not sub-pixel rendering is used during data deposition.
        periodic: Whether the pixelization will span the domain boundaries.
        POutput: Whether or not it should print runtime info, on by default as defined by init().
    Returns:
        An FRB
    """
    
    if isinstance(bounds, str):
        if bounds.upper() == "XY":
            bounds = (source.ds.domain_left_edge[0], source.ds.domain_right_edge[0], source.ds.domain_left_edge[1], source.ds.domain_right_edge[1])
        elif bounds.upper() == "YZ":
            if source.ds.dimensionality == 2:
                raise TypeError("Cannot use Z axis for bounds on 2D data")
            bounds = (source.ds.domain_left_edge[1], source.ds.domain_right_edge[1], source.ds.domain_left_edge[2], source.ds.domain_right_edge[2])
        elif bounds.upper() == "XZ":
            if source.ds.dimensionality == 2:
                raise TypeError("Cannot use Z axis for bounds on 2D data")
            bounds = (source.ds.domain_left_edge[0], source.ds.domain_right_edge[0], source.ds.domain_left_edge[2], source.ds.domain_right_edge[2])
        else:
            raise TypeError("Unknown fullDomain argument")
    
    if level > source.ds.max_level:
        if ((POutput == None) and (PLog == True)) or (POutput == True): print("Level provided is too high, defaulting to maximum possible level")
        level = source.ds.max_level
    
    resLevel =  (source.ds.domain_dimensions[0]*2**(level), source.ds.domain_dimensions[1]*2**(level))

    if len(bounds) != 4:
        raise TypeError("Bounds bust be a length 4 array of form (Xmin, Xmax, Ymin, Ymax)")
    frb = yt.FixedResolutionBuffer(source, bounds, (resLevel[1],resLevel[0]),antialias, periodic)
    if ((POutput == None) and (PLog == True)) or (POutput == True): print("Generating fixed resolution buffer")
    return frb

#Goes through making a slice and turning it into an FRB directly and automatically    
def genSliceArray(ds,axis,coords,bounds,field, resLevel=None,sourceS=None,center=None,fieldParameters=None, antialias=True, YTArray=True, POutput=None, periodic=False):
    """
    A function that automatically slices, makes an FRB and then an array
    
    Parameters:
        ds: The dataset to use
        axis: The axis to cut across. If the given axis is not possible, it will raise an error. Can be 0,1,2s for x,y,z. Cylinders can only be cut along 2
        coordinate: The coordinates where the slice will be made
        bounds: Bounds are the min and max in the image plane that we want our image to cover.  It's in the order of (xmin, xmax, ymin, ymax), where the coordinates are all in the appropriate code units. Can automatically use the maximum by providing strins "XY", "YZ" and "XZ"
        field: The field to use if making an array
        center: The center to pass to fields that require it
        sourceS: Source to use instead of all of ds
        fieldParameters: A dictionary of field parameters than can be accessed by derived fields.
        resLevel: The resolution level to use to find the size of the image
        antialias: Determines whether or not sub-pixel rendering is used during data deposition.
        YTArray: Whether The function returns a numpy array or a YT array (with units). It will return a YT array by default
        POutput: Whether or not it should print runtime info, on by default as defined by init().
    Returns:
        An array, either a numpy array or a YT array
    """
    if resLevel == None:
        resLevel = ds.max_level
    
    if isinstance(bounds, str):
        if bounds == "full domain":
            if axis == 0:
                bounds = (ds.domain_left_edge[2], ds.domain_right_edge[2], ds.domain_left_edge[1], ds.domain_right_edge[1])
            elif axis == 1:
                bounds = (ds.domain_left_edge[2], ds.domain_right_edge[2], ds.domain_left_edge[0], ds.domain_right_edge[0])
            elif axis == 2:
                bounds = (ds.domain_left_edge[0], ds.domain_right_edge[0], ds.domain_left_edge[1], ds.domain_right_edge[1])
            else:
                raise ValueError("Invalid axis")
        else:
            raise TypeError("Unknown bounds command, must be full domain")
                
    slc = sliceData(ds, axis, coords, center, fieldParameters, source=sourceS)
    frb = makeFRB(slc, bounds, resLevel,antialias,periodic=periodic)
    
    if YTArray == True:
        if ((POutput == None) and (PLog == True)) or (POutput == True): print( "Returning YTArray")
        frbf = frb[field]
        if axis == 1:
            frbf = np.rot90(frbf, 1)
        return frbf#I think this part doesn't work. Double check how YT makes arrays vs numpy arrays
        
    else:
        if ((POutput == None) and (PLog == True)) or (POutput == True): print( "Returning Numpy array")
        npArray = np.array(frb[field]) #Something here doesn't seem quite right
        if axis == 1:
            npArray = np.rot90(npArray, 1)
        return npArray              

#Generates a line out array of the provided data set    
#add some prints
#check whether the points provided meet the dimensionality required
def genLineOutArray(ds,axis,points,field,fieldParameters=None,dataSource=None,YTArray=True, POutput=None):
    """
    A function that creates a lineout and returns it as an array. Works on 1D,2D and 3D
    
    Parameters:
        ds: The dataset to lineout
        axis: The axis to slice along. If the given axis is not possible, it will raise an error. Can be 0,1,2 for x,y,z
        points: tuple of floats. The (plane_x, plane_y) coordinates at which to cast the ray. Note that this is in the plane coordinates: so if you are casting along x, this will be (y, z).  If you are casting along y, this will be (z, x).  If you are casting along z, this will be (x, y).
        fieldParameters: A dictionary of field parameters than can be accessed by derived fields.
        dataSource: Draw the selection from the provided data source rather than all data associated with the data_set
        YTArray: Whether this function returns a numpy array or a YT array (with units). It will return a YT array by default
        POutput: Whether or not it should print runtime info, on by default as defined by init().
    Returns:
        An array, either a numpy array or a YT array
    """
    if ((ds.dimensionality == 1) or (ds.dimensionality == 2)) and (axis == 2):
        dimensionalityForError = str(ds.dimensionality) + "D"
        raise TypeError("Cannot use Z-axis for %s geometry" % dimensionalityForError)
    if(ds.dimensionality == 1) and (axis == 1):
        raise TypeError("Cannot use the Y-axis for 1D geometry")
    if checkPossibleLocation(ds,points) == False:
        raise TypeError("Points specified are out of bounds")
    lno = ds.ortho_ray(axis, points,field_parameters=fieldParameters,data_source=dataSource)
    if YTArray == True:
        return lno[field]
    else:
        lno = np.array(lno[field])
        return lno
        
#Extract a region of raw data from 3D FLASH yt dataset.
#dims en cm; center en cm, YTarray of dimension 3.        
def getRegion(ds, level = 0, dims=None, center = None,gZones=None,POutput=None):
    """
    A function that extracts a region of data from a 3D Flash dataset.
    
    Parameters:
        ds: The dataset from which to extract.
        level: The refinement level to be used
        dims: Dimensions of the area (in cm)
        center: Center of the area
        POutput: Whether or not it should print runtime info, on by default as defined by init().
    Returns:
        A covering grid
    """
    if center is None: 
        raise TypeError("Center need to be indicated")
    elif not isinstance(center, yt.units.yt_array.YTArray): # Verify input center_yt is a YT Array, e.g. not just a NumPy array
        raise TypeError("Input 'center_yt' must be a YT Array (of type yt.units.yt_array.YTArray). E.g. try center_yt = ds.arr([0,0,3],'cm')")
    if dims is None: 
        raise TypeError("Dimensions need to be indicated (cm)")
    elif not isinstance(center, yt.units.yt_array.YTArray): # Verify input center_yt is a YT Array, e.g. not just a NumPy array
        raise TypeError("Input 'dims' must be a YT Array (of type yt.units.yt_array.YTArray).")
    # Extract covering-grid cube with given center and edge length
    
    level = checkPossibleLevel(ds,level)
    
    
    le = center - dims/2.0 # YT left edge for cube
    re = center + dims/2.0 # YT right edge for cube
    if ((POutput == None) and (PLog == True)) or (POutput == True): 
        print( "Left edge of cube: ", le)
        print( "Right edge of cube: ", re)
    cellwid0 = ds.domain_width/ds.domain_dimensions # Cell width in YT units, at level=0 refinement (lrefine=1)
    if ((POutput == None) and (PLog == True)) or (POutput == True): print( "Cell width at level 0 refinement: ", cellwid0)
    dimensions = np.array((re - le)/cellwid0) * 2**level + 1 # Desired array dimensions ([int, int, int]) of covering grid depends on cell resolution
    if ((POutput == None) and (PLog == True)) or (POutput == True): print("Desired array dimensions: ", dimensions)
    if gZones != None:
        cg = ds.covering_grid(level=level, left_edge=le, dims=dimensions, num_ghost_zones=gZones) # Covering grid cube
    else:
        cg = ds.covering_grid(level=level, left_edge=le, dims=dimensions)

    return cg

#fine cell size at the given refinement level
def cellSize(ds,refinement=0,POutput=None):
    """
    A function that finds cell size at a given refinement level
    
    Parameters:
        ds: The dataset from which to extract.
        refinement: The refinement level to use when calculating size.
        POutput: Whether or not it should print runtime info, on by default as defined by init().
    Returns:
        A YT array
    """
    #check refinement doesn't go over ds.max_level
    #if so, print problem and default to maximum
    
    tref = checkPossibleLevel(ds,refinement)
    
    cellWid = ds.domain_width/(ds.domain_dimensions*2**tref)
    if ds.dimensionality == 1:
        if ((POutput == None) and (PLog == True)) or (POutput == True): print("Cell width at refinement level ", tref, "  is ", cellWid[0])
        return cellWid[0]
    elif ds.dimensionality == 2:
        if ((POutput == None) and (PLog == True)) or (POutput == True): print("Cell width at refinement level ", tref, "  is ", (cellWid[0],cellWid[1]))
        return (cellWid[0],cellWid[1])
    if ((POutput == None) and (PLog == True)) or (POutput == True): print("Cell width at refinement level ", tref, "  is ", cellWid)
    return cellWid



#check the center +- dims/2 is within bounds    
def cubeDim(ds,center,dims,level,POutput=None):
    """
    A function that finds the dimension of a cube for 3D data
    
    Parameters:
        ds: The dataset from which to extract.
        center: Center of the cube (in cm)
        dims: Dimension of the cube (in cm)
        level: Level of refinement to use, cannot go above the maximum number of units
        unit: Unit to use with the YT array for center and dims (REMOVE)
        POutput: Whether or not it should print runtime info, on by default as defined by init().
    Returns:
        A YT array
    """
    if len(dims) != 3:
        raise TypeError("Bounds must have 3 arguments")
    if len(center) != 3:
        raise TypeError("Center must have 3 arguments")
    if checkPossibleLocation(ds,center,POutput) == False:
        raise TypeError("Center is in an invalid location")
    if ds.dimensionality != 3:
        dimensionalityForError = str(ds.dimensionality) + "D"
        raise TypeError("Cannot cube data with a dimensionality of %s" % dimensionalityForError)
    Ccenter = ds.arr(center,"cm")
    Cdims = ds.arr(dims,"cm")
    le = Ccenter - Cdims/2.0 # YT left edge for cube
    re = Ccenter + Cdims/2.0 # YT right edge for cube
    if ((POutput == None) and (PLog == True)) or (POutput == True): 
        print("Left edge of cube: ", le)
        print("Right edge of cube: ", re)
    cellWid = cellSize(ds)
    if ((POutput == None) and (PLog == True)) or (POutput == True): print("Cell width at refinement level ", tref, "  is ", cellWid)
    dimensions = np.array((re - le)/cellWid) * 2**level + 1 # Desired array dimensions ([int, int, int]) of covering grid depends on cell resolution
    return dimensions
    
#add smoothed vs unsmooted option    
#add checker for above maximum refinement
#add tonumpyarray function
def genCubeCoveringGrid(ds, center, dims, level,field=None,gZones=0,pbar=True,fieldParameters=None, smoothed=False,POutput=None):
    """
    A function that produces an array of the given 3D data
    
    Parameters:
        ds: The dataset from which to extract.
        center: Center of the cube
        dims: Number of cells along each axis
        level: Level of resolution to use
        unit: Unit to use with the YT array
        le: Left edge of dataset
        field: List of fields
        gZones: Number of padding ghost zones
        pbar: Whether to use pbar or not. True by default
        fieldParameters: A dictionary of field parameters than can be accessed by derived fields.
        smoothed: Whether the grid should be smoothed or not, off by default.
        POutput: Whether or not it should print runtime info, on by default as defined by init().
    Returns:
        A covering grid
    """
    Ccenter = ds.arr(center,"cm")
    Cdims = ds.arr(dims,"cm")
    le = Ccenter - Cdims/2.0 # YT left edge for cube
    if len(le) != 3:
        raise TypeError("le must have 3 arguments")
    if checkPossibleLocation(ds,le, POutput) == False:
        print("Left edge defined is out of bounds, defaulting to left edge")
        le = ds.domain_left_edge
    dimension = cubeDim(ds, center, dims, level, unit)
    if ((POutput == None) and (PLog == True)) or (POutput == True): print("Dimension provided by cubeDim: ", dimension)
    ogDim = ds.domain_dimensions
    cellWid = cellSize(ds)
    i = 0
    """
    while(i < 3):
        dimension[i] = ogDim[i]
        i = i+1
        print dimension
            
    dimension = dimension - np.array(2*cellWid)
    """
    if smoothed == True:
      cg = ds.smoothed_covering_grid(level=level, left_edge=le, dims=dimension, fields=field, num_ghost_zones=gZones,use_pbar=pbar,field_parameters=fieldParameters)
      return cg
    cg = ds.covering_grid(level=level, left_edge=le, dims=dimension, fields=field, num_ghost_zones=gZones,use_pbar=pbar,field_parameters=fieldParameters)
    return cg
    
def genProj(ds, field, axis, weightField=None, center=None, dataSource=None, method="integrate", fieldParameters=None,maxLevel=None,POutput=None):
    """
    A function that just generates the projection
    Parameters:
        ds: The dataset from which to extract.
        field: Field to use
        axis: The Axis to cut across. If the given axis is not possible, it will raise an error. Can be 0,1,2 for x,y,z
        weightField: Multiplier for the field 
        center: Center for the fields that need it
        fieldParameters: Values to be passed as field parameters that can be accessed by generated fields.
        dataSource: Specific source for data
        method: The method of projection (integrate,mip,sum)
        maxLevel: The maximum allowed level
        POutput: Whether or not it should print runtime info, on by default as defined by init().
    Returns:
        Projection object which can be accessed or tuned into an FRB
        
"""
    if (ds.geometry.upper() == "CYLINDRICAL") and (axis != 2):
        raise TypeError("You may only project cylinders along axis 2")    
    if ((ds.dimensionality == 2) and (axis == 2)) and (ds.geometry.upper() != "CYLINDRICAL"):
        dimensionalityForError = str(ds.dimensionality) + "D " + ds.geometry
        raise TypeError("Cannot use the Z-axis for %s geometry" % dimensionalityForError)
    elif (ds.dimensionality == 1):
        raise TypeError("Cannot project 1D geometry")
    if ((POutput == None) and (PLog == True)) or (POutput == True): print("Generating projection along axis", axis)
    prj = ds.proj(field,axis,weight_field=weightField,center=center,data_source=dataSource,method=method,field_parameters=fieldParameters,max_level=maxLevel)    
    
    
    
def genProjFRB(ds,field,axis,width,resolution,weightField=None,center=None, fieldParameters=None,dataSource=None,maxLevel=None,centerFRB=None,height=None,periodic=False,YTArray=True,method='integrate',POutput=None):
    """
    A function that generates an FRB based on a projection. It's a little slow though
    
    Parameters:
        ds: The dataset from which to extract.
        field: Field to use
        axis: The Axis to cut across. If the given axis is not possible, it will raise an error. Can be 0,1,2 for x,y,z
        width: Width of the data
        resolution: resolution for the FRB
        weightField: Multiplier for the field 
        center: Center for the fields that need it
        fieldParameters:
        dataSource: Specific source for data
        maxLevel: The maximum allowed level
        centerFRB: Center for the FRB
        height: Height of the data
        periodic: Whether the pixelization will span the domain boundaries.
        YTArray: Whether the array should be numpy or YT. YT by default
        method: The method of projection (integrate,mip,sum)
`       POutput: Whether or not it should print runtime info, on by default as defined by init().
    Returns:
        An array, either numpy or YT
"""
    if (ds.geometry.upper() == "CYLINDRICAL") and (axis != 2):
        raise TypeError("You may only project cylinders along axis 2")    
    if ((ds.dimensionality == 2) and (axis == 2)) and (ds.geometry.upper() != "CYLINDRICAL"):
        dimensionalityForError = str(ds.dimensionality) + "D " + ds.geometry
        raise TypeError("Cannot use the Z-axis for %s geometry" % dimensionalityForError)
    elif (ds.dimensionality == 1):
        raise TypeError("Cannot project 1D geometry")
    prj = ds.proj(field,axis,weight_field=weightField,center=center,data_source=dataSource,method=method,field_parameters=fieldParameters,max_level=maxLevel)
    frb = prj.to_frb(width, resolution, centerFRB, height, periodic)
    if ((POutput == None) and (PLog == True)) or (POutput == True): print("Generating projection along %s axis" % axis)
    if YTArray == False:
        if ((POutput == None) and (PLog == True)) or (POutput == True): print("Returning Numpy array style")
        npArray = toNumpyArray(frb, field)
        return npArray
    else:
        if ((POutput == None) and (PLog == True)) or (POutput == True): print("Returning YT style FRB")
        return frb[field]



def fetchField(source, field, loc='flash', YTArray=True):
    """
    A function that pulls the data from a field
    
    Parameters:
        source: The source to pull data from (what is the source? An YT object that is not ds)
        field: The field to pull from
        loc: The location of the field, it's "Flash" by default
        YTArray: Whether the data should be returned as a YT array or a numpy Array, YT array by default
    Returns:
        An array, either YT or numpy
    """
    if YTArray == True:
        return source[(loc,field)]
    else:
        return np.array(source[(loc,field)])
    

def singlePoint(ds,p,fieldParameters=None,dataSource=None,POutput=None,field=None):
    """
    A function that pulls the data from a single provided point
    
    Parameters:
        ds: The dataset to pull the point from
        p: the point to pull from as a list
        fieldParameters: A dictionary of field parameters than can be accessed by derived fields.
        dataSource: Specific source for data
        POutput: Whether or not it should print runtime info, on by default as defined by init().
    Returns:
        A point object
    """
    if ds.dimensionality == 1:
        if isinstance(p,float) or isinstance(p, int):
            p = (p,0.0,0.0)
    if ds.dimensionality == 2:
        if len(p) != 3:
            p = p.append(0.0)
    if checkPossibleLocation(ds,p,POutput) == False:
        raise TypeError("Coordinates out of bounds")
    pnt = ds.point(p,field_parameters=fieldParameters,data_source=dataSource)
    if field != None:
        return pnt[field]
    return pnt

#add documentation, used mainly internally for other things    


##############INTERNAL FUNCTION USED BY OTHER FUNCTIONS##################
def checkPossibleLocation(ds,r, POutput=None):
        """
        A function used internally that checks whether a provided location is within the bounds of domain.
        
        Parameters:
            ds: ds to check the bounds for.
            r: The location to check. Form (x,y,z)
            POutput: Whether or not it should print runtime info, on by default as defined by init().
        Returns:
            True or False depending on whether the location is within bounds
        
        
        """
        le = ds.domain_left_edge
        re = ds.domain_right_edge
        if len(r) == 3:
            if (le[0] <= r[0]) and (re[0] >= r[0]):
                if (le[1] <= r[1]) and (re[1] >= r[1]):
                    if (le[2] <= r[2]) and (re[2] >= r[2]):
                        if ((POutput == None) and (PLog == True)) or (POutput == True): print("All values in bounds")
                        return True   
                    if ((POutput == None) and (PLog == True)) or (POutput == True): print("Z value out of bounds")    
                    return False
                if ((POutput == None) and (PLog == True)) or (POutput == True): print("Y value out of bounds") 
                return False
            if ((POutput == None) and (PLog == True)) or (POutput == True): print("X value out of bounds")
            return False  
        elif len(r) == 2:
            if (le[0] <= r[0]) and (re[0] >= r[0]):
                if (le[1] <= r[1]) and (re[1] >= r[1]):
                    if ((POutput == None) and (PLog == True)) or (POutput == True): print("All values in bounds")
                    return True
                if ((POutput == None) and (PLog == True)) or (POutput == True): print("Y value out of bounds") 
                return False
            if ((POutput == None) and (PLog == True)) or (POutput == True): print("X value out of bounds")
            return False
            
def basic2DPlot(frb,axisLabel=None,name=None,colormap="default",title=None,POutput=None,interpolation=None,extent=None):
    """
    A function for constructing 1D plots. Saves the image to the directory currently being worked in
    
    Parameters:
        frb: The dataset to graph. It can either take a YT-style FRB or a numpy array.
        colormap: The colormap to use for the plot. It can accept the following colormaps:
            AUTUMN
            BONE
            COOL
            FLAG
            GRAY
            HOT
            HSV
            INFERNO
            JET
            MAGMA
            NIPY SPECTRAL
            PINK
            PLASMA
            PRISM
            SPRING
            SUMMER
            VIRIDIS
            WINTER
            
        axis: Words to put on the X and Y axis as well as on the colorbar. Pass it a last with the order ["Xaxis","YAxis","Colorbar"]
        name: The name to save the file as. It will by default save as a .png if not specified.
        title: Title at the top of the plot.
        interpolation: The interpolation method to use, supports none, nearest, bilinear, bicubic, spline16, splice36,hanning, hamming, hermite, kaiser, quadric, catrom, gaussian, bessel, mitchell, sinc and lanczos
        POutput: Whether or not it should print runtime info, on by default as defined by init().
    
    """
    
    if isinstance(frb,yt.ImageArray):
        if ((POutput == None) and (PLog == True)) or (POutput == True): print("Loading up data as YT FRB")
        show = np.array(frb)
        plt.figure(figsize=(12,8))
        plt.imshow(show, interpolation=interpolation, extent=extent)
    elif isinstance(frb,np.ndarray): 
        if ((POutput == None) and (PLog == True)) or (POutput == True): print("Loading up data as numpy array")
        plt.figure(figsize=(12,8))
        plt.imshow(frb, interpolation=interpolation, extent=extent)
    else:
        raise TypeError("Unknown array type, please pass either a numpy array or a YT array")
    
    
    colormapSet(colormap)    
    q = plt.colorbar()
    if ((POutput == None) and (PLog == True)) or (POutput == True):
        print("Setting color map to ",  colormap.upper())
    if axisLabel != None:
        if ((POutput == None) and (PLog == True)) or (POutput == True): print("Setting axes" )
        yax = plt.ylabel(axisLabel[0])
        plt.xlabel(axisLabel[1])
        if len(axisLabel) >= 3:
            q.set_label(axisLabel[2])          
        yax.set_rotation(0)
    if title != None:
        if ((POutput == None) and (PLog == True)) or (POutput == True): print("Setting title")
        plt.suptitle(title)
    #figure out how to change the colormap to be logarithmic
    if name != None:
        plt.savefig(name)
    #plt.close()
    
def basic1DPlot(frb, format="-b", axis=None, name=None, title=None,XKCD=False,POutput=None,scaleX=None,scaleY=None,axisLimX=None, axisLimY=None):
    """
    A function for constructing 1D plots. Saves the image to the directory currently being worked in
    
    Parameters:
        frb: The dataset to graph. It can either take a YT-style FRB or a numpy array.
        format: The color and appearance of the lines, as well as any point markers. Defaults to solid blue lines with no markers. It can take the following:
            Line types:
                '-'	    solid line style
                '--'	    dashed line style
                '-.'	    dash-dot line style
                ':'	    dotted line style
            Line colors:
                'b'	    blue
                'g'	    green
                'r'	    red
                'c'	    cyan
                'm'	    magenta
                'y'	    yellow
                'k'	    black
                'w'	    white
            Point markers:
                '.'	    point marker
                ','	    pixel marker
                'o'	    circle marker
                'v'	    triangle_down marker
                '^'	    triangle_up marker
                '<'	    triangle_left marker
                '>'	    triangle_right marker
                '1'	    tri_down marker
                '2'	    tri_up marker
                '3'	    tri_left marker
                '4'	    tri_right marker
                's'	    square marker
                'p'	    pentagon marker
                '*'	    star marker
                'h'	    hexagon1 marker
                'H'	    hexagon2 marker
                '+'	    plus marker
                'x'	    x marker
                'D'	    diamond marker
                'd'	    thin_diamond marker
                '|'	    vline marker
                '_'	    hline marker
                
        axis: Words to put on the X and Y axis. Pass it a last with the order ["Xaxis","YAxis"]
        name: The name to save the file as. It will by default save as a .png if not specified.
        title: Title at the top of the plot.
        XCKD: Turn on XCKD style drawing. Included just for laughs
        POutput: Whether or not it should print runtime info, on by default as defined by init().
    
    
    """
    if XKCD == True:
        if ((POutput == None) and (PLog == True)) or (POutput == True): print("Setting XKCD draw mode")
        plt.xkcd()
        
    if axisLimX != None:
        if ((POutput == None) and (PLog == True)) or (POutput == True): print("Setting scale for X-axis")
        plt.xlim(axisLimX)
        
    if axisLimY != None:
        if ((POutput == None) and (PLog == True)) or (POutput == True): print("Setting scale for Y-axis")
        plt.ylim(axisLimY)
        
    if isinstance(frb,yt.ImageArray):
        if ((POutput == None) and (PLog == True)) or (POutput == True): print("Loading up data as YT FRB")
        plt.figure()
        plt.plot(np.array(frb),format)   
    else:
        if ((POutput == None) and (PLog == True)) or (POutput == True): print("Loading up data as numpy array")
        plt.figure()
        plt.plot(frb,format)
        
    if axis != None:
        if ((POutput == None) and (PLog == True)) or (POutput == True): print("Setting axes")
        plt.ylabel(axis[0])
        plt.xlabel(axis[1])
    if title != None:
        if ((POutput == None) and (PLog == True)) or (POutput == True): print("Setting title")
        pltsuptitle(title)    
    if name != None:
        plt.savefig(name)
    #plt.close()
    
    
##############INTERNAL FUNCTION USED BY OTHER FUNCTIONS##################
def colormapSet(colormap, POutput=None):
    """
    A function used internally that sets the colormap of 2D plots.
    
    Parameters:
        colormap: The name of the colormap to use, It defaults to the default one within the function using it.
    Returns:
        None
 
    """
    if colormap.upper() == "DEFAULT":           
        pass
    elif colormap.upper() == "AUTUMN":          
        plt.autumn()
    elif colormap.upper() == "BONE":            
        plt.bone()
    elif colormap.upper() == "COOL":            
        plt.cool()       
    elif colormap.upper() == "COPPER":          
        plt.copper()
    elif colormap.upper() == "FLAG":            
        plt.flag()
    elif (colormap.upper() == "GRAY") or (colormap.upper() == "GREY"):            
        plt.gray()
    elif colormap.upper() == "HOT":             
        plt.hot()
    elif colormap.upper() == "HSV":             
        plt.hsv()
    elif colormap.upper() == "INFERNO":         
        plt.inferno()
    elif colormap.upper() == "JET":             
        plt.jet()
    elif colormap.upper() == "MAGMA":           
        plt.magma()
    elif colormap.upper() == "NIPY SPECTRAL":   
        plt.nipy_spectral()
    elif colormap.upper() == "PINK":            
        plt.pink()
    elif colormap.upper() == "PLASMA":          
        plt.plasma()
    elif colormap.upper() == "PRISM":           
        plt.prism()
    elif colormap.upper() == "SPRING":          
        plt.spring()
    elif colormap.upper() == "SUMMER":          
        plt.summer()
    elif colormap.upper() == "VIRIDIS":         
        plt.viridis()
    elif colormap.upper() == "WINTER":          
        plt.winter()
    else:                                       
        print("Invalid colormap, using default")
 
 
def getFullDomainExtent(ds,axis,POutput=None):
    if axis == 0:      
        q = [ds.domain_left_edge[2],ds.domain_right_edge[2],ds.domain_left_edge[1],ds.domain_right_edge[1]]
        if ((POutput == None) and (PLog == True)) or (POutput == True): print("Full domain for 0 is: ", q)
    elif axis == 1:
        q = [ds.domain_left_edge[0],ds.domain_right_edge[0],ds.domain_left_edge[2],ds.domain_right_edge[2]]
        if ((POutput == None) and (PLog == True)) or (POutput == True): print("Full domain for 1 is: ", q)
    elif axis == 2:
        q = [ds.domain_left_edge[0],ds.domain_right_edge[0],ds.domain_left_edge[1],ds.domain_right_edge[1]]
        if ((POutput == None) and (PLog == True)) or (POutput == True): print("Full domain for 2 is: ", q)
    else:
        raise ValueError('Invalid axis to slice along')
    return q
 
##############INTERNAL FUNCTION USED BY OTHER FUNCTIONS##################
def checkPossibleLevel(ds,rflevel,POutput=None):
    """
    A function used internally that checks whether a provided refinement level is possible
    
    Parameters:
        ds: ds to check the max level for.
        rflevel: The refinement level to check against the max.
        POutput: Whether or not it should print runtime info, on by default as defined by init().
    Returns:
        Either the provided refinement level if it is possible, or returns the maximum possible refinement if the provided level is too high.
 
    """
    if ds.max_level < rflevel:
        if ((POutput == None) and (PLog == True)) or (POutput == True):print("Provided refinement level is too high. Defaulting to the maximum", ds.max_level)
        return ds.max_level
    else:
        if ((POutput == None) and (PLog == True)) or (POutput == True):
            print("Provided refinement level within bounds")
        return rflevel

def plot2DDS(ds, axis, coords, bounds, field, logScale=False,resLevel=None, sourceS=None, center=None, fieldParameters=None, antialias=True, POutput=None, periodic=False, axisLabels=None, name=None, colormap="default", title=None, interpolation=None,ext=None):
    print("This may take a few moments")        
    data = genSliceArray(ds,axis,coords,bounds,field,resLevel,sourceS,center,fieldParameters,antialias,POutput=POutput,periodic=periodic,YTArray=False)
    if logScale:
        data =np.log10(data)
    if ext == None:
        ext = getFullDomainExtent(ds,axis,POutput)
    if axisLabels == None:
        if axis == 0:
            ax1 = "Y"
            ax2 = "Z"
        elif axis == 1:
            ax1 = "X"
            ax2 = "Z"
        elif axis == 2:
            ax1 = "X"
            ax2 = "Y"
        r1 = str(ax1+ ": "+ str(ds.unit_system['length']))
        r2 = str(ax2+ ": "+ str(ds.unit_system['length']))
        axisLabels = [r1,r2]
    if title == None:
        if field == 'dens': title = str("Density: " + getUnits(ds,field))
        elif field == 'temp': title = str("Temperature: "+ getUnits(ds,field))
        elif field == 'tele': title = str("Electron Temperature: "+ getUnits(ds,field))
        elif field == 'trad': title = str("Radiation Temperature: "+ getUnits(ds,field))
        elif field == 'velx': title = str("X Velocity: "+ getUnits(ds,field))
        elif field == 'vely': title = str("Y Velocity: "+ getUnits(ds,field))
        elif field == 'velz': title = str("Z Velocity: "+ getUnits(ds,field))
        elif field == 'tion': title = str("Ion Temperature: "+ getUnits(ds,field))
        elif field == 'magx': title = str(":"+ getUnits(ds,field))
        elif field == 'magy': title = str(":"+ getUnits(ds,field))
        elif field == 'magz': title = str(":"+ getUnits(ds,field))
        elif field == 'magp': title = str(":"+ getUnits(ds,field))
        elif field == 'pres': title = str("Pressure:"+ getUnits(ds,field))
    
    return (data,axisLabels,name,colormap,title,interpolation,ext)

def plot1DDS(ds,axis,points,field, logScale=False,fieldParameters=None,dataSource=None,POutput=None,format='-b',axisNames=None,name=None,title=None,XKCD=False,axisLimX=None, axisLimY=None):
    print("This may take a moment")
    data = genLineOutArray(ds,axis,points,field,fieldParameters,dataSource,YTArray=False,POutput=POutput)   
    if logScale:
        data =np.log10(data)

    if axisLimX == 'full domain':
       axisLimX = singleAxisDomain(ds,axis)
    
    basic1DPlot(data, format=format, axis=axisNames, name=name, title=title,XKCD=XKCD,POutput=POutput,axisLimX=axisLimX, axisLimY=axisLimY)

def singleAxisDomain(ds,axis,POutput=None):    
    if isinstance(axis, int):
        if ((POutput == None) and (PLog == True)) or (POutput == True):print("Finding domain limits for axis", axis)
        return (ds.domain_left_edge[axis],ds.domain_right_edge[axis])
    else:
        raise TypeError("Axis input can be 0,1 or 2")

def getUnits(ds,field,POutput=None):
    if ((POutput == None) and (PLog == True)) or (POutput == True): print("Fetching units for field")
    temp = str(ds.unit_system['temperature'])
    length = str(ds.unit_system['length'])
    time = str(ds.unit_system['time'])
    mass = str(ds.unit_system['mass'])
    
    if field == 'dens': return str(mass+"/"+length+"^3")
    elif field == 'temp': return temp
    elif field == 'tele': return temp
    elif field == 'trad': return temp
    elif field == 'velx': return str(length+'/'+time)
    elif field == 'vely': return str(length+'/'+time)
    elif field == 'velz': return str(length+'/'+time)
    elif field == 'tion': return temp
    elif field == 'magx': return str('sqrt('+mass+")/"+'(sqrt('+length+")^"+time+")")
    elif field == 'magy': return str('sqrt('+mass+")/"+'(sqrt('+length+")^"+time+")")
    elif field == 'magz': return str('sqrt('+mass+")/"+'(sqrt('+length+")^"+time+")")
    elif field == 'magp': return str(mass+"/("+length+"*"+time+"^2)")
    elif field == 'pres': return str(mass+"/("+length+"*"+time+"^2)")
    else: 
        print("Unknown field, either I haven't added it or you're an idiot who can't type")
        pass
    

#TODO:
#Change around unit display on plot, have it place the equastion in the title next to the field, as well as displaying the axis it is (X,Y,Z) and the unit of length
#Do it for 1D plot too
#Add in the streamplot stuff that is in gorgon



#We should probably add more catches for 1d-ten-t errors, IBM errors and layer 8 issues ID:10T

