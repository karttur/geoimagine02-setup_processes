'''
Created on 20 Jan 2021

@author: thomasgumbricht
'''

from ancillary import ProcessAncillary

from modis import ProcessModis

from postgresdb import ManageProcess, ManageRegion, ManageMODIS

from gis import kt_gis as ktgis

from os import path, makedirs

from sys import exit

from params import JsonParams

from base64 import b64encode

import netrc

from setup_processes import PGsession

from postgresdb import ManageAncillary

from params.layers import VectorLayer
 

def DbConnect(db):
    '''
    '''
    # the HOST must exist in the .netrc file in the users home directory
    HOST = 'karttur'
    
    # Retrieve login and password from the .netrc file
    secrets = netrc.netrc()
    
    # Authenticate username, account and password 
    username, account, password = secrets.authenticators( HOST )
    
    # Encode the password before sending it
    password = b64encode(password.encode())

    # Create a query dictionary for connecting to the Postgres server  
    query = {'db':db, 'user':username, 'pswd':password}
    
    return query

class ProcessProcess:  
    """"class for processes defining other processes"""  
    
    def __init__(self, pp): 
        """"The constructor requires an instance of the main process, and the json object defining the process to setup"""  
                
        self.pp = pp
 
        self.session = ManageProcess(self.pp.postgresdb.db)
                
        if self.pp.process.processid == 'addrootproc':
            
            queryD = {'rootprocid':self.pp.process.parameters.rootprocid, 
                    'title':self.pp.process.parameters.title,
                  'label':self.pp.process.parameters.label,
                  'creator':self.pp.userproject.userid}

            self.session._ManageRootProcess(self.pp.process, queryD)
        
        elif self.pp.process.processid == 'addsubproc':
            
            queryD = {'rootprocid':self.pp.process.parameters.rootprocid,
                      'subprocid':self.pp.process.parameters.subprocid, 
                      'title':self.pp.process.parameters.title,
                      'label':self.pp.process.parameters.label,
                      'version':self.pp.process.parameters.version,
                      'minuserstratum':self.pp.process.parameters.minuserstratum,
                      'creator':self.pp.userproject.userid}
        
            self.session._ManageSubProcess(self.pp.process, queryD)
        
        else:
            
            exitstr = 'subprocess %s not defined in manageprocess' %(self.process.processid)
            
            exit( exitstr )

class ProcessDefaultRegions:
    '''
    '''  
    def __init__(self, pp ):
        
        self.pp = pp
 
        self.session = ManageRegion(self.pp.postgresdb.db)
        
        #direct to subprocess
        if self.pp.process.processid == 'RegionCategories':
                        
            self.session._InsertRegionCat(self.pp.process)
            
        elif self.pp.process.processid == 'DefaultRegionFromCoords':
        
            self._DefaultRegionFromCoords()
            
        elif self.pp.process.processid == 'DefaultRegionFromVector':
            
            self._DefaultRegFromVec()
            
        elif self.process.proc.processid == 'linkregionswrs':
            if self.process.proj.projectid == "karttur" and self.process.proj.system == 'system': 
                self.LinkAllRegionsToWRS()
        elif self.process.proc.processid == 'linkregionsmodtiles':
            if self.process.proj.projectid == "karttur" and self.process.proj.system == 'system': 
                self.LinkAllRegionsToMODIS()
            else:
                exit('Only superuser can link wrs to regions')
        else:
            exitstr = 'No process %s under Processregion' %(self.process.processid)
            exit(exitstr)
            
    def _DefaultRegionRegister(self,layer):
        '''
        '''
        
        # Get the projection
        projection = ktgis.GetVectorProjection(layer.FPN)
        
        #Set lonlat projection
        lonlatproj = ktgis.MjProj()
        
        lonlatproj.SetFromEPSG(4326)
        
        # Get the boundary
        boundsD = ktgis.GetFeatureBounds(layer.FPN,'REGIONID')  
        
        if len(boundsD) != 1:
        
            exitstr = 'Default regions must consist on only one (1) feature (polygon or multipolygon): %s' %(layer.FPN)
            
            exit(exitstr)                   
                
        k = list(boundsD)[0]
                
        layer._SetBounds(projection.epsg,boundsD[k][0], boundsD[k][1], boundsD[k][2], boundsD[k][3] )  
                
        #Get the corners in lonlat
        llD = ktgis.ReprojectBounds(layer.BoundsPtL,projection.proj_cs,lonlatproj.proj_cs)

        queryD = {'regionid': self.pp.process.parameters.regionid, 
                  'regionname': self.pp.process.parameters.regionname,
                  'parentid': self.pp.process.parameters.parentid, 
                  'regioncat': self.pp.process.parameters.regioncat, 
                  'parentcat': self.pp.process.parameters.parentcat,
                  'title': self.pp.process.parameters.title, 
                  'label':self.pp.process.parameters.label,
                  'epsg':self.pp.procsys.dstepsg}
        
        if ' ' in queryD['regionid'] or ' ' in queryD['parentid']:
            
            exit('regionid or parentid with whuite space in setup_process_process') 
            
        self.session._InsertDefRegion(layer, queryD, boundsD[self.pp.process.parameters.regionid], llD, self.pp.process.overwrite, self.pp.process.delete )
         
    def _DefaultRegionFromCoords(self):
        '''
        '''

        for locus in self.pp.dstLayerD:
            
            for datum in self.pp.dstLayerD[locus]:
                
                for comp in self.pp.dstLayerD[locus][datum]:
                                        
                    layer = self.pp.dstLayerD[locus][datum][comp]

                    #The destination region must be forced,this is because the locus to be created did not exists when checking for the default locus
 
                    layer.locus.locus = self.pp.process.parameters.regionid.lower()
                      
                    layer.locus.path = self.pp.process.parameters.regionid.lower()
                                                                      
                    layer._SetPath()
                    
                    fieldDD = self._SetfieldD()
                       
                    layer.CreateVectorAttributeDef(fieldDD)
                    
                    layer._SetBounds(self.pp.procsys.dstepsg, 
                                     self.pp.process.parameters.minx, 
                                     self.pp.process.parameters.miny, 
                                     self.pp.process.parameters.maxx, 
                                     self.pp.process.parameters.maxy)

                    projection = ktgis.MjProj()
                    
                    projection.SetFromEPSG(self.pp.procsys.dstepsg)

                    if not layer._Exists() or self.pp.process.overwrite:
                        
                        ktgis.CreateESRIPolygonPtL(layer.FPN, layer.fieldDefL, layer.BoundsPtL, projection.proj_cs, self.pp.process.parameters.regionid)          
                    
                    self._DefaultRegionRegister(layer)
                    
                    
    def _DefaultRegFromVec(self):
        '''
        '''
        
        # dstLayerD and srcLayerD are almost identical
        for locus in self.pp.dstLayerD:
            
            for datum in self.pp.dstLayerD[locus]:
                
                for comp in self.pp.dstLayerD[locus][datum]:
                                                                
                    srcLayer = self.pp.srcLayerD[locus][datum][comp]
                    
                    if not path.isfile(srcLayer.FPN):
                        
                        exitstr = 'No source layer in _DefaultRegFromVec', srcLayer.FPN
                        
                        exit(exitstr)
                        
                    p = self.pp.process.parameters
                    
                    fieldL = [p.vector_db_id, p.vector_db_name,
                               p.vector_db_category, p.vector_db_parentid, 
                               p.vector_db_parentcat, p.vector_db_stratum,
                               p.vector_db_title, p.vector_db_label]
                    
                    fieldD = ktgis.GetFeatureAttributeList(srcLayer.FPN, fieldL, p.vector_db_id)
                    
                    if not fieldD:
                        exit('setup_process_class: fieldD failed in _DefaultRegFromVec')
                    
                    for key in fieldD:
                        
                        # Convert the field data to a dict 
                        params = ['regionid', 'regionname', 'regioncat', 'stratum', 'parentid', 'parentcat', 'title', 'label']
                        
                        values = [ str(fieldD[key][p.vector_db_id]).lower().replace(' ', '-'), 
                                    str(fieldD[key][p.vector_db_name]), 
                                    str(fieldD[key][p.vector_db_category].lower()), 
                                    int(fieldD[key][p.vector_db_stratum]), 
                                    str(fieldD[key][p.vector_db_parentid]).lower().replace(' ', '-'),
                                    str(fieldD[key][p.vector_db_parentcat].lower()),
                                    str(fieldD[key][p.vector_db_title]),
                                    str(fieldD[key][p.vector_db_label]) ]
                        
                        d = dict(zip(params, values))
                        
                        # Replace the process class parameter with the dict      
                        self.pp.process.parameters = lambda:None
                        
                        for k,v in d.items():
                            
                            setattr(self.pp.process.parameters, k, v) 
                            
                        fieldDD = self._SetfieldD()
                                                                    
                        regionid = self.pp.process.parameters.regionid
                        
                        #Construct the locus for this region
                        locusD = {'locus':regionid,'path':regionid}
                        
                        # Recreate the composition
                        compDstCopy = self.pp.dstLayerD[locus][datum][comp].comp
                        
                        # Set layerid and prefix to "defreg"
                        compDstCopy.layerid = compDstCopy.prefix = 'defreg'
                        
                        # Set content ot roi (region of interest)
                        compDstCopy.content = 'roi'
                        
                        # Reset the compid
                        compDstCopy._SetCompid()
                        
                        # Recreate the vector Layer
                        dstLayer = VectorLayer(compDstCopy, locusD, self.pp.dstPeriod.datumD[datum])
                    
                        dstLayer.CreateVectorAttributeDef(fieldDD)
                        
                        fieldname = p.vector_db_id
                    
                        valueLL = [[fieldD[key][p.vector_db_id]]]
                        
                        if not dstLayer._Exists() or self.pp.process.overwrite: #or overwrite
                           
                            ktgis.ExtractFeaturesToNewDS(srcLayer.FPN, dstLayer.FPN, fieldname,valueLL, dstLayer.fieldDefL)
                         
                        self._DefaultRegionRegister(dstLayer)
                            
                        '''
                            fieldname = 'REGIONID'
                            
                            #Get the epsg and bounds
                            boundsD = ktgis.GetFeatureBounds(dstLayer.FPN,fieldname) 
                            
                            if len(boundsD) != 1:
                            
                                exitstr = 'Default regions must consist on only one (1) feature (polygon or multipolygon): %s' %(dstLayer.FPN)
                                
                                exit(exitstr)
                            
                            projection = ktgis.GetVectorProjection(dstLayer.FPN)
                            
                            k = list(boundsD)[0]
                            
                            bounds = boundsD[k]
            
                            dstLayer._SetBounds(projection.epsg,boundsD[k][0], boundsD[k][1], boundsD[k][2], boundsD[k][3] )  
                  
                            _DefaultRegionRegister(self,dstLayer, projection)
                            
                            #Set lonlat projection
                            
                            lonlatproj = ktgis.MjProj()
                            
                            lonlatproj.SetFromEPSG(4326)
                            
                            #Get the corners in lonlat
                            
                            llD = ktgis.ReprojectBounds(dstLayer.BoundsPtL,projection.proj_cs,lonlatproj.proj_cs)
            
            
                            title = label = 'default region %s' %(regionid)
                            
                            query = {'regionname':regionname,'regioncat':regioncat, 'parentid':parentid, 'parentcat':parentcat,'regionid':regionid, 'title':title,'label':label,'epsg':projection.epsg}
                            
                            session._InsertDefRegion(self.process, dstLayer, query, bounds, llD )
                        '''                          

    
    def _SetfieldD(self):
        ''' Set the fields for default region layers
        '''
        
        fieldDD = {}
        
        fieldDD['REGIONID'] = {'name':'REGIONID', 'type':'string','width':32,
                               'precision':0,'transfer':'constant','source':self.pp.process.parameters.regionid }
        
        fieldDD['NAME'] = {'name':'NAME', 'type':'string','width':64,
                           'precision':0,'transfer':'constant','source':self.pp.process.parameters.regionname }
        
        fieldDD['CATEGORY'] = {'name':'CATEGORY', 'type':'string','width':32,
                               'precision':0,'transfer':'constant','source':self.pp.process.parameters.regioncat }
        
        fieldDD['STRATUM'] = {'name':'STRATUM', 'type':'integer','width':4,
                              'precision':0,'transfer':'constant','source':self.pp.process.parameters.stratum }
        
        fieldDD['PARENTID'] = {'name':'PARENTID', 'type':'string','width':32,
                               'precision':0,'transfer':'constant','source':self.pp.process.parameters.parentid }
        
        fieldDD['PARENTCAT'] = {'name':'PARENTCAT', 'type':'string','width':32,
                                'precision':0,'transfer':'constant','source':self.pp.process.parameters.parentcat }
                
        return fieldDD
    
def SetupProcessesRegions(docpath, projFN, db):
    '''
    Setup processes
    '''
    
    srcFP = path.join(path.dirname(__file__),docpath)
    
    projFPN = path.join(srcFP,projFN)
    
    # Get the full path to the project text file
    dirPath = path.split(projFPN)[0]
    
    if not path.exists(projFPN):
        
        exitstr = 'EXITING, project file missing: %s' %(projFPN)
        
        exit( exitstr )
    
    infostr = 'Processing %s' %(projFPN)
    
    print (infostr)
    
    # Open and read the text file linking to all json files defining the project
    with open(projFPN) as f:
        
        jsonL = f.readlines()
    
    # Clean the list of json objects from comments and whithespace etc    
    jsonL = [path.join(dirPath,x.strip())  for x in jsonL if len(x) > 10 and x[0] != '#']
        
    # Get the user and password for connecting to the db
    query = DbConnect(db)

    # Connect to the Postgres Server
    session = PGsession(query)
        
    ProcPar = JsonParams(session)
    
    processL = []
    
    #Loop over all json files and create Schemas and Tables
    for jsonObj in jsonL:
                
        processL.append( ProcPar._JsonObj(jsonObj) )
        
    # Close the db connection for getting processes and user
    session._Close()
    
    for processD in processL:

        for k in range(len(processD)):
            
            print ('    ',k, processD[k]) 

            if processD[k]['PP'].rootprocid == 'manageprocess':
    
                ProcessProcess(processD[k]['PP'])
                    
            elif processD[k]['PP'].rootprocid == 'ManageRegion':
    
                #ProcessDefaultRegions(db, process, self.procsys, self.userproject, self.userid, self.usercat, self.stratum)
                ProcessDefaultRegions(processD[k]['PP'])
                
            elif processD[k]['PP'].rootprocid == 'Ancillary':
    
                session = ManageAncillary(db)
                
                ProcessAncillary(processD[k]['PP'], session)
                
                session._Close()
                
            elif processD[k]['PP'].rootprocid == 'MODISProc':
    
                session = ManageMODIS(db)
                
                ProcessModis(processD[k]['PP'], session)
                
                session._Close()
                
            else:
                
                print (processD[k]['PP'].rootprocid)
                print (processD[k]['PP'].subprocid)

def ModisTileCoords(db, verbose = 1):
    ''' Create the MODIS defaut tiling system
    '''
            
    #Open the db session for MODIS
    session = ManageMODIS(db)
    
    SINproj = ktgis.MjProj()
    
    SINproj.SetFromProj4('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs')
    
    LatLonproj = ktgis.MjProj()
    
    LatLonproj.SetFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')
    
    ptL = []
    
    for lon in range(360):
        
        ptL.append((lon-180,90))
        
    for lat in range(180):
        
        ptL.append((180,-1*(lat-90)))
        
    for lon in range(360):
        
        ptL.append((-1*(lon-180),-90))
        
    for lat in range(180):
        
        ptL.append((-180,lat-90))
        

    worldgeom = ktgis.ShapelyPolyGeom(ptL)
    
    worldgeom.ShapelyToOgrGeom()
    
    worldgeom.GeoTransform(LatLonproj,SINproj)
    
    worldgeom.OgrGeomToShapely()
    
    home = path.expanduser("~")
    
    tarShpFP =  path.join(path.dirname(__file__),'data')

    if not path.exists(tarShpFP):
        
        makedirs(tarShpFP)
        
    FN = 'modtiles-multi_karttur_global_epsg6842.shp'
    
    tarShpFPN = path.join(tarShpFP,FN) 
    
    fieldDefD = {'type':'string','transfer':'constant','source':'globe','width':8}
    
    fieldDefL = [ktgis.FieldDef('name',fieldDefD)]
    
    ktgis.CreateESRIPolygonGeom(tarShpFPN, fieldDefL, worldgeom, SINproj.proj_cs, 'globe')
    
    # Create a shape file for all individual tiles in SIN proj
    
    FN = 'modtiles-single_karttur_global_epsg6842.shp'
    
    tarShpFPN = path.join(tarShpFP,FN)
    
    tarDS,tarLayer = ktgis.ESRICreateDSLayer(tarShpFPN, SINproj.proj_cs, 'polygon', 'tiles', fieldDefL)
    
    # Create a shape file for all individual tiles in Geographic coordinates
    FN = 'modtiles_karttur_global_0.shp'
    
    tarShpFPN = path.join(tarShpFP,FN)
    
    tarDSLonLat,tarLayerLonLat = ktgis.ESRICreateDSLayer(tarShpFPN, LatLonproj.proj_cs, 'polygon', 'tiles', fieldDefL)
    
    #create a region with all tiles
    tlen = 20015109.3539999984204769
    
    tlen /= 18

    for h in range(36):
        
        minx = tlen*(18-36)+h*tlen
        
        maxx = minx+tlen
        
        for v in range(18):
            
            maxy = tlen*(9-18)+(18-v)*tlen
            
            miny = maxy-tlen
            
            ptL = [(minx,maxy),(maxx,maxy),(maxx,miny),(minx,miny)]
            
            tilegeom = ktgis.ShapelyMultiPointGeom(ptL)
            
            #convert to ogr
            tilegeom.ShapelyToOgrGeom()
            
            #write target feature 
            tilegeom.GeoTransform(SINproj,LatLonproj)
            
            tilegeom.OgrGeomToShapely()
            
            coordL = []
            
            for point in [ptgeom for ptgeom in tilegeom.shapelyGeom]:
                
                coordL.extend([list(point.coords)[0][0],list(point.coords)[0][1]])
                
            ullon, ullat, urlon, urlat, lrlon, lrlat, lllon, lllat = coordL  
              
            tilepoly = ktgis.ShapelyPolyGeom([(minx, maxy), (maxx, maxy), (maxx, miny), (minx,miny)])
            
            #Test if this tile is inside the globe
            if tilepoly.shapelyGeom.intersects(worldgeom.shapelyGeom): 
                
                if h < 10:
                    
                    htile = 'h0%s' %(h)
                    
                else:
                    
                    htile = 'h%s' %(h)
                    
                if v < 10:
                    
                    vtile = 'v0%s' %(v)
                    
                else:
                    
                    vtile = 'v%s' %(v)
                    
                hvtile = '%s%s' %(htile,vtile)
                
                polytilegeom = ktgis.ShapelyPolyGeom(ptL)
                
                polytilegeom.ShapelyToOgrGeom()
                
                fieldDefD = {'type':'string','transfer':'constant','source':hvtile,'width':8}
                
                fieldDefL = [ktgis.FieldDef('name',fieldDefD)]
                
                #create target feature
                tarFeat = ktgis.ogrFeature(tarLayer)
                                
                tarFeat.CreateOgrFeature(polytilegeom.ogrGeom, fieldDefL) 
                
                if h == 17:
                    
                    pass
                
                else:
                    
                    #to be correct 5 points are needed and also the lat must be fitted
                    if h < 18 and ullon > 0:
                        
                        ullon = -180
                        
                    if h < 18 and lllon > 0:
                        
                        lllon = -180

                    if h < 18 and urlon > 0:
                        
                        urlon = -180

                    if h < 18 and lrlon > 0:
                        
                        lrlon = -180

                    if h > 18 and urlon < 0:
                        
                        urlon = 180
 
                    if h > 18 and lrlon < 0:
                        
                        lrlon = 180

                    if h > 18 and ullon < 0:
                        
                        ullon = 180
 
                    if h > 18 and lllon < 0:
                        
                        lllon = 180
                        
                    if hvtile == 'h24v01':
                        
                        urlon = 180
                        
                    if hvtile == 'h24v16':
                        
                        lrlon = 180  
                        
                    if hvtile == 'h11v01':
                        
                        ullon = -180
                        
                    if hvtile == 'h11v16':
                        
                        lllon = -180 
                        
                if ullon > urlon:
                    
                    print ('ERROR','ullon > urlon',hvtile,ullon,urlon)

                if lllon > lrlon:
                    
                    print ('ERROR','lllon > lrlon',hvtile, lllon, lrlon)

                #
                polytilegeom = ktgis.ShapelyPolyGeom([(ullon, ullat), (urlon, urlat), (lrlon, lrlat), (lllon,lllat)])
                
                polytilegeom.ShapelyToOgrGeom()
                
                #polytilegeom.GeoTransform(SINproj,LatLonproj)
                
                #create target feature
                tarLonLatFeat = ktgis.ogrFeature(tarLayerLonLat)
                
                tarLonLatFeat.CreateOgrFeature(polytilegeom.ogrGeom, fieldDefL)
                
                west,south,east,north = polytilegeom.shapelyGeom.bounds

                session._InsertModisTileCoord(hvtile,h,v,minx,maxy,maxx,miny,west,south,east,north,ullat,ullon,lrlon,lrlat,urlon,urlat,lllon,lllat)

                query = {'system':'system','table':'regions','h':h,'v':v,'hvtile':hvtile,'regionid':'global','regioncat':'global','regiontype':'default','delete':False}
                
                session._InsertModisRegionTile(query)
                
    tarDS.CloseDS()
    
    tarDSLonLat.CloseDS()
    
    session._Close()
    
    if verbose > 1:
        
        print ('Check the shaoe file',tarShpFPN)
    
    return (tarShpFPN)


def Ease2NTileCoords(db, verbose = 1):
    ''' Create the MODIS defaut tiling system
    '''
            
    #Open the db session for Ease2
    #session = ManageEase2(db)
    
    Ease2Nproj = ktgis.MjProj()
    
    Ease2Nproj.SetFromProj4('+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
    
    LatLonproj = ktgis.MjProj()
    
    LatLonproj.SetFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')
    
    ptL = []
    
    for lon in range(360):
        
        ptL.append((lon-180,0))
        
    for lat in range(90):
        
        ptL.append((180,lat))
        
    for lon in range(360):
        
        ptL.append((-1*(lon-180),0))
        
    for lat in range(90):
        
        ptL.append((-180,lat))
        

    worldgeom = ktgis.ShapelyPolyGeom(ptL)
    
    worldgeom.ShapelyToOgrGeom()
    
    worldgeom.GeoTransform(LatLonproj,Ease2Nproj)
    
    worldgeom.OgrGeomToShapely()
    
    home = path.expanduser("~")
    
    tarShpFP =  path.join(path.dirname(__file__),'data')

    if not path.exists(tarShpFP):
        
        makedirs(tarShpFP)
        
    FN = 'ease2ntiles-multi_karttur_northpole_epsg6931.shp'
    
    tarShpFPN = path.join(tarShpFP,FN) 
    
    fieldDefD = {'type':'string','transfer':'constant','source':'globe','width':8}
    
    fieldDefL = [ktgis.FieldDef('name',fieldDefD)]
    
    ktgis.CreateESRIPolygonGeom(tarShpFPN, fieldDefL, worldgeom, Ease2Nproj.proj_cs, 'northpole')
    
    # Create a shape file for all individual tiles in SIN proj
    
    FN = 'ease2ntiles-single_karttur_northpole_epsg6931.shp'
    
    tarShpFPN = path.join(tarShpFP,FN)
    
    tarDS,tarLayer = ktgis.ESRICreateDSLayer(tarShpFPN, Ease2Nproj.proj_cs, 'polygon', 'tiles', fieldDefL)
    
    # Create a shape file for all individual tiles in Geographic coordinates
    FN = 'ease2ntiles_karttur_northpole_0.shp'
    
    tarShpFPN = path.join(tarShpFP,FN)
    
    tarDSLonLat,tarLayerLonLat = ktgis.ESRICreateDSLayer(tarShpFPN, LatLonproj.proj_cs, 'polygon', 'tiles', fieldDefL)
    
    #create a region with all tiles
    tlen = 20015109.3539999984204769
    
    tlen /= 18

    for h in range(36):
        
        minx = tlen*(18-36)+h*tlen
        
        maxx = minx+tlen
        
        for v in range(18):
            
            maxy = tlen*(9-18)+(18-v)*tlen
            
            miny = maxy-tlen
            
            ptL = [(minx,maxy),(maxx,maxy),(maxx,miny),(minx,miny)]
            
            tilegeom = ktgis.ShapelyMultiPointGeom(ptL)
            
            #convert to ogr
            tilegeom.ShapelyToOgrGeom()
            
            #write target feature 
            tilegeom.GeoTransform(Ease2Nproj,LatLonproj)
            
            tilegeom.OgrGeomToShapely()
            
            coordL = []
            
            for point in [ptgeom for ptgeom in tilegeom.shapelyGeom]:
                
                coordL.extend([list(point.coords)[0][0],list(point.coords)[0][1]])
                
            ullon, ullat, urlon, urlat, lrlon, lrlat, lllon, lllat = coordL  
              
            tilepoly = ktgis.ShapelyPolyGeom([(minx, maxy), (maxx, maxy), (maxx, miny), (minx,miny)])
            
            #Test if this tile is inside the globe
            if tilepoly.shapelyGeom.intersects(worldgeom.shapelyGeom): 
                
                if h < 10:
                    
                    htile = 'h0%s' %(h)
                    
                else:
                    
                    htile = 'h%s' %(h)
                    
                if v < 10:
                    
                    vtile = 'v0%s' %(v)
                    
                else:
                    
                    vtile = 'v%s' %(v)
                    
                hvtile = '%s%s' %(htile,vtile)
                
                polytilegeom = ktgis.ShapelyPolyGeom(ptL)
                
                polytilegeom.ShapelyToOgrGeom()
                
                fieldDefD = {'type':'string','transfer':'constant','source':hvtile,'width':8}
                
                fieldDefL = [ktgis.FieldDef('name',fieldDefD)]
                
                #create target feature
                tarFeat = ktgis.ogrFeature(tarLayer)
                                
                tarFeat.CreateOgrFeature(polytilegeom.ogrGeom, fieldDefL) 
                
                if h == 17:
                    
                    pass
                
                else:
                    
                    #to be correct 5 points are needed and also the lat must be fitted
                    if h < 18 and ullon > 0:
                        
                        ullon = -180
                        
                    if h < 18 and lllon > 0:
                        
                        lllon = -180

                    if h < 18 and urlon > 0:
                        
                        urlon = -180

                    if h < 18 and lrlon > 0:
                        
                        lrlon = -180

                    if h > 18 and urlon < 0:
                        
                        urlon = 180
 
                    if h > 18 and lrlon < 0:
                        
                        lrlon = 180

                    if h > 18 and ullon < 0:
                        
                        ullon = 180
 
                    if h > 18 and lllon < 0:
                        
                        lllon = 180
                              
                if ullon > urlon:
                    
                    print ('ERROR','ullon > urlon',hvtile,ullon,urlon)

                if lllon > lrlon:
                    
                    print ('ERROR','lllon > lrlon',hvtile, lllon, lrlon)

                #
                polytilegeom = ktgis.ShapelyPolyGeom([(ullon, ullat), (urlon, urlat), (lrlon, lrlat), (lllon,lllat)])
                
                polytilegeom.ShapelyToOgrGeom()
                
                #polytilegeom.GeoTransform(SINproj,LatLonproj)
                
                #create target feature
                tarLonLatFeat = ktgis.ogrFeature(tarLayerLonLat)
                
                tarLonLatFeat.CreateOgrFeature(polytilegeom.ogrGeom, fieldDefL)
                
                west,south,east,north = polytilegeom.shapelyGeom.bounds

                #session._InsertModisTileCoord(hvtile,h,v,minx,maxy,maxx,miny,west,south,east,north,ullat,ullon,lrlon,lrlat,urlon,urlat,lllon,lllat)

                query = {'system':'system','table':'regions','h':h,'v':v,'hvtile':hvtile,'regionid':'global','regioncat':'global','regiontype':'default','delete':False}
                
                #session._InsertModisRegionTile(query)
                
    tarDS.CloseDS()
    
    tarDSLonLat.CloseDS()
    
    #session._Close()
    
    if verbose > 1:
        
        print ('Check the shaoe file',tarShpFPN)
    
    return (tarShpFPN)

def Ease2NTileCoords6391(db, verbose = 1):
    ''' Create the MODIS defaut tiling system
    '''
            
    #Open the db session for Ease2
    #session = ManageEase2(db)
    
    Ease2Nproj = ktgis.MjProj()
    
    Ease2Nproj.SetFromProj4('+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
    
    LatLonproj = ktgis.MjProj()
    
    LatLonproj.SetFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')
    
    ptL = []
    
    minx = -15725.34707351
    miny = 9009951.03827368
    maxx = 974.37914072
    maxy = 558277.55721523
    
    centerx = (minx+maxx)/2
    
    centery = (miny+maxy)/2
    
    
    print ('centerx', centerx)
    
    print ('centery', centery)
    
    BAlle
    
    '''
    -15725.34707351 9009951.03827368
    974.37914072 558277.55721523
'''
    
    for lon in range(360):
        
        ptL.append((lon-180,0))
        
    for lat in range(90):
        
        ptL.append((180,lat))
        
    for lon in range(360):
        
        ptL.append((-1*(lon-180),0))
        
    for lat in range(90):
        
        ptL.append((-180,lat))
        

    worldgeom = ktgis.ShapelyPolyGeom(ptL)
    
    worldgeom.ShapelyToOgrGeom()
    
    worldgeom.GeoTransform(LatLonproj,Ease2Nproj)
    
    worldgeom.OgrGeomToShapely()
    
    home = path.expanduser("~")
    
    tarShpFP =  path.join(path.dirname(__file__),'data')

    if not path.exists(tarShpFP):
        
        makedirs(tarShpFP)
        
    FN = 'ease2ntiles-multi_karttur_northpole_epsg6931.shp'
    
    tarShpFPN = path.join(tarShpFP,FN) 
    
    fieldDefD = {'type':'string','transfer':'constant','source':'globe','width':8}
    
    fieldDefL = [ktgis.FieldDef('name',fieldDefD)]
    
    ktgis.CreateESRIPolygonGeom(tarShpFPN, fieldDefL, worldgeom, Ease2Nproj.proj_cs, 'northpole')
    
    # Create a shape file for all individual tiles in SIN proj
    
    FN = 'ease2ntiles-single_karttur_northpole_epsg6931.shp'
    
    tarShpFPN = path.join(tarShpFP,FN)
    
    tarDS,tarLayer = ktgis.ESRICreateDSLayer(tarShpFPN, Ease2Nproj.proj_cs, 'polygon', 'tiles', fieldDefL)
    
    # Create a shape file for all individual tiles in Geographic coordinates
    FN = 'ease2ntiles_karttur_northpole_0.shp'
    
    tarShpFPN = path.join(tarShpFP,FN)
    
    tarDSLonLat,tarLayerLonLat = ktgis.ESRICreateDSLayer(tarShpFPN, LatLonproj.proj_cs, 'polygon', 'tiles', fieldDefL)
    
    #create a region with all tiles
    tlen = 20015109.3539999984204769
    
    tlen /= 18

    for h in range(36):
        
        minx = tlen*(18-36)+h*tlen
        
        maxx = minx+tlen
        
        for v in range(18):
            
            maxy = tlen*(9-18)+(18-v)*tlen
            
            miny = maxy-tlen
            
            ptL = [(minx,maxy),(maxx,maxy),(maxx,miny),(minx,miny)]
            
            tilegeom = ktgis.ShapelyMultiPointGeom(ptL)
            
            #convert to ogr
            tilegeom.ShapelyToOgrGeom()
            
            #write target feature 
            tilegeom.GeoTransform(Ease2Nproj,LatLonproj)
            
            tilegeom.OgrGeomToShapely()
            
            coordL = []
            
            for point in [ptgeom for ptgeom in tilegeom.shapelyGeom]:
                
                coordL.extend([list(point.coords)[0][0],list(point.coords)[0][1]])
                
            ullon, ullat, urlon, urlat, lrlon, lrlat, lllon, lllat = coordL  
              
            tilepoly = ktgis.ShapelyPolyGeom([(minx, maxy), (maxx, maxy), (maxx, miny), (minx,miny)])
            
            #Test if this tile is inside the globe
            if tilepoly.shapelyGeom.intersects(worldgeom.shapelyGeom): 
                
                if h < 10:
                    
                    htile = 'h0%s' %(h)
                    
                else:
                    
                    htile = 'h%s' %(h)
                    
                if v < 10:
                    
                    vtile = 'v0%s' %(v)
                    
                else:
                    
                    vtile = 'v%s' %(v)
                    
                hvtile = '%s%s' %(htile,vtile)
                
                polytilegeom = ktgis.ShapelyPolyGeom(ptL)
                
                polytilegeom.ShapelyToOgrGeom()
                
                fieldDefD = {'type':'string','transfer':'constant','source':hvtile,'width':8}
                
                fieldDefL = [ktgis.FieldDef('name',fieldDefD)]
                
                #create target feature
                tarFeat = ktgis.ogrFeature(tarLayer)
                                
                tarFeat.CreateOgrFeature(polytilegeom.ogrGeom, fieldDefL) 
                
                if h == 17:
                    
                    pass
                
                else:
                    
                    #to be correct 5 points are needed and also the lat must be fitted
                    if h < 18 and ullon > 0:
                        
                        ullon = -180
                        
                    if h < 18 and lllon > 0:
                        
                        lllon = -180

                    if h < 18 and urlon > 0:
                        
                        urlon = -180

                    if h < 18 and lrlon > 0:
                        
                        lrlon = -180

                    if h > 18 and urlon < 0:
                        
                        urlon = 180
 
                    if h > 18 and lrlon < 0:
                        
                        lrlon = 180

                    if h > 18 and ullon < 0:
                        
                        ullon = 180
 
                    if h > 18 and lllon < 0:
                        
                        lllon = 180
                              
                if ullon > urlon:
                    
                    print ('ERROR','ullon > urlon',hvtile,ullon,urlon)

                if lllon > lrlon:
                    
                    print ('ERROR','lllon > lrlon',hvtile, lllon, lrlon)

                #
                polytilegeom = ktgis.ShapelyPolyGeom([(ullon, ullat), (urlon, urlat), (lrlon, lrlat), (lllon,lllat)])
                
                polytilegeom.ShapelyToOgrGeom()
                
                #polytilegeom.GeoTransform(SINproj,LatLonproj)
                
                #create target feature
                tarLonLatFeat = ktgis.ogrFeature(tarLayerLonLat)
                
                tarLonLatFeat.CreateOgrFeature(polytilegeom.ogrGeom, fieldDefL)
                
                west,south,east,north = polytilegeom.shapelyGeom.bounds

                #session._InsertModisTileCoord(hvtile,h,v,minx,maxy,maxx,miny,west,south,east,north,ullat,ullon,lrlon,lrlat,urlon,urlat,lllon,lllat)

                query = {'system':'system','table':'regions','h':h,'v':v,'hvtile':hvtile,'regionid':'global','regioncat':'global','regiontype':'default','delete':False}
                
                #session._InsertModisRegionTile(query)
                
    tarDS.CloseDS()
    
    tarDSLonLat.CloseDS()
    
    #session._Close()
    
    if verbose > 1:
        
        print ('Check the shaoe file',tarShpFPN)
    
    return (tarShpFPN)