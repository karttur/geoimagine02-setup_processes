'''
Created on 21 feb. 2018
Updated on 4 jan. 2021

@author: thomasgumbricht
'''
import psycopg2

from params import JsonParams, Struct

from base64 import b64decode

from sys import exit

#from os import path
from pprint import pprint

from postgresdb import ManageProcess, ManageLayout, ManageRegion

class ProcessProcess:  
    """"class for processes defining other processes"""  
    
    def __init__(self, db, process, userid, usercat, stratum): 
        """"The constructor requires an instance of the main process, and the json object defining the process to setup"""  
                
        self.session = ManageProcess(db)
        
        if process.processid == 'addrootproc':
            
            queryD = {'rootprocid':process.parameters.rootprocid, 
                    'title':process.parameters.title,
                  'label':process.parameters.rootprocid,
                  'creator':userid}

            self.session._ManageRootProcess(process, queryD)
        
        elif process.processid == 'addsubproc':
            
            queryD = {'rootprocid':process.parameters.rootprocid,
                      'subprocid':process.parameters.subprocid, 
                    'title':process.parameters.title,
                  'label':process.parameters.rootprocid,
                  'version':process.parameters.version,
                  'minuserstratum':process.parameters.minuserstratum,
                  'creator':userid}
        
            self.session._ManageSubProcess(process, queryD)
        
        else:
            
            exitstr = 'subprocess %s not defined in manageprocess' %(self.process.processid)
            
            exit( exitstr )
            
class PGsession:
    """Connect to postgres server"""  
     
    def __init__(self, query):
        """Connect to selected database"""
        
        query['pswd'] = b64decode(query['pswd']).decode('ascii')
        
        conn_string = "host='localhost' dbname='%(db)s' user='%(user)s' password='%(pswd)s'" %query
                
        self.conn = psycopg2.connect(conn_string)
        
        self.cursor = self.conn.cursor()
            
    def _DictToSelect(self, queryD):
        '''
        Converts a dictionary to Select statement 
        '''
        selectL = []
        for key in queryD:
            #statement = key operator value
            statement = ' %(col)s %(op)s \'%(val)s\'' %{'col':key.replace('#',''), 'op':queryD[key]['op'], 'val':queryD[key]['val']}
            selectL.append(statement)
        self.select_query = "WHERE %(where)s" %{'where':' AND '.join(selectL)}  
        return self.select_query
    
    def _SelectRootProcess(self,queryD):
        '''
        '''
        self.cursor.execute("SELECT rootprocid, minuserstratum FROM process.subprocesses WHERE subprocid = '%(subprocid)s';" %queryD)
        
        record = self.cursor.fetchone()
        
        return record
    
    def _SelectUserCred(self, queryD):
        '''
        '''
                
        sql = "SELECT userid, usercat, stratum FROM userlocale.users WHERE userid = '%(user)s';" %queryD
        
        self.cursor.execute(sql)    
        
        self.userid, self.usercat, self.stratum = self.cursor.fetchone()
                
    def _SelectTractDefRegion(self, queryD):
        '''
        '''
        #First check if this region is itself a defregion
        
        sql = "SELECT regionid FROM regions.defregions WHERE regionid = '%(tract)s';" %queryD
        
        self.cursor.execute(sql)
        
        rec = self.cursor.fetchone()
        
        if rec != None:
        
            return (rec[0], 'D')

        sql = "SELECT parentid FROM regions.tracts WHERE tractid = '%(tract)s';" %queryD
        
        self.cursor.execute(sql) 
        
        rec = self.cursor.fetchone()
        
        if rec == None:
            
            return rec
        
        return (rec[0], 'T')
        
    def _SelectProcessSystem(self, queryD, paramL):
        ''' Select system for this process
        '''
        
        queryD['cols'] = " ,".join(paramL)
        
        sql = "SELECT %(cols)s FROM process.procsys WHERE subprocid = '%(subprocid)s' and system = '%(system)s';" %queryD

        self.cursor.execute(sql)    
        
        record = self.cursor.fetchone()
        
        if record == None:
                        
            self.cursor.execute("SELECT srcsystem, dstsystem, srcdivision, dstdivisio FROM process.procsys WHERE subprocid = '%(subprocid)s' and system = '*';" %queryD)    

            record = self.cursor.fetchone() 
        
        if record == None:
            
            SNULLE
        
        return dict(zip(paramL,record)) 
    
    def _SingleSearch(self,queryD, paramL,  table, schema, pq = False):
        #self._GetTableKeys(schema, table)
        selectQuery = {}
        for item in queryD:

            if isinstance(queryD[item],dict):
                #preset operator and value
                selectQuery[item] = queryD[item]
            else:
                selectQuery[item] = {'op':'=', 'val':queryD[item]}
        wherestatement = self._DictToSelect(selectQuery)  
        cols =  ','.join(paramL)
        selectQuery = {'schema':schema, 'table':table, 'select': wherestatement, 'cols':cols}
        
        query = "SELECT %(cols)s FROM %(schema)s.%(table)s %(select)s" %selectQuery
        if pq:
            print ('SingleSearch query',query)
        self.cursor.execute(query)
        self.records = self.cursor.fetchone()
        return self.records
          
    def _SelectMulti(self,queryD, paramL, table, schema):  
        ''' Select multiple records from any schema.table
        '''

        selectQuery = {}
        for item in queryD:

            if isinstance(queryD[item],dict):
                #preset operator and value
                selectQuery[item] = queryD[item]
            else:
                selectQuery[item] = {'op':'=', 'val':queryD[item]}
        wherestatement = self._DictToSelect(selectQuery) 

        if len(paramL) == 1:
            cols = paramL[0]
        else:
            cols =  ','.join(paramL)
        selectQuery = {'schema':schema, 'table':table, 'select': wherestatement, 'cols':cols}      
        query = "SELECT %(cols)s FROM %(schema)s.%(table)s %(select)s" %selectQuery
 
        self.cursor.execute(query)
        
        self.records = self.cursor.fetchall()
        
        return self.records    
            
    def ReadRunJson(self, jsonObj, db):
        ''' Read and run json object
        '''
        
        PP = Params(jsonObj)
        
        self.params = params.params
        
        self.jsonParams = params.jsonParams
        
        # Convert dict {userproject} to Struct
        self.userproject = Struct(self.jsonParams['userproject'])
        
        if (self.params.process[0].verbose):
            
            pprint(self.jsonParams)
            
        # Get the processes as listed in the jsonObject
        for process in self.params.process:
                        
            self.verbose = process.verbose
            
            self.delete = process.delete
            
            self.overwrite = process.overwrite
            
            queryD ={'subprocid':process.processid}
            
            
            # Get the root process id for the process
            #processid = self._SingleSearch(queryD, paramL, system, 'layers', True)
            
            record = self._SelectRootProcess(queryD)
            
            if record == None:
                
                exitstr = 'Exiting - The process %s is not registered in the Framework db' %(process.processid)
            
                exit(exitstr)
                
            self.rootprocid, self.minuserstratum = record
            
            if self.verbose > 1:
                
                print ('process parameters')
                
                #pprint (vars(process))
                
                print (process.parameters.rootprocid[0])
                
            # Get the system associated with this process
                
            queryD['system'] = self.userproject.system
            
            paramL = ('srcsystem', 'dstsystem', 'srcdivision', 'dstdivision', 'srcepsg', 'dstepsg')
            
            record = self._SelectProcessSystem(queryD, paramL)
            
            if record == None:
                
                SNULLEBULLE
                exit()
              
            # Convert the system setting to a class object called system
            self.procsys = lambda: None

            for k, v in record.items():
                
                setattr(self.procsys, k, v)
                               
            if self.rootprocid == 'manageprocess':

                ProcessProcess(db, process, self.userid, self.usercat, self.stratum)
                
            elif self.rootprocid == 'ManageRegion':

                ProcessDefaultRegions(db, process, self.procsys, self.userproject, self.userid, self.usercat, self.stratum)
                
                
                
            else:
                
                print (self.rootprocid)
                SNYLTARE
                
    def _Close(self):
        
        self.cursor.close()
        
        self.conn.close()
        
class ProcessDefaultRegionsOld:
    '''
    '''  
    def __init__(self, db, process, procsys, userproject, userid, usercat, stratum):

        self.session = ManageRegion(db)
        
        self.process = process
        
        self.procsys = procsys
        
        self.userproject = userproject
        
        #direct to subprocess
        if self.process.processid == 'RegionCategories':
            
            self.session._InsertRegionCat(self.process)
            
        #elif self.process.proc.processid == 'defaultregion':
        elif self.process.processid == 'DefaultRegionFromCoords':
        
            self._DefaultRegion(self.process, self.procsys, self.userproject)
            
        elif self.process.proc.processid == 'defaultregionfromvector':
            self._DefaultRegFromVec(session)
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
                    
    def _DefaultRegion(self,process, procsys, userproject):
        '''
        '''
        
        locus = 'hej'
        for locus in self.process.dstLayerD:
            
            for datum in self.process.dstLayerD[locus]:
                
                for comp in self.process.dstLayerD[locus][datum]:
                    
                    query = self.process.proc.paramsD
                    fieldDD = self._SetfieldD( query['regionid'], query['regionname'], query['regioncat'], query['stratum'], query['parentid'], query['parentcat'])
                    layer = self.process.dstLayerD[locus][datum][comp]

                    #The destination region must be forced,this is because the locus to be created did not exists when chekcing for the feault locus
                    if self.verbose:
                        print ('        forcing region from to',locus, self.process.proc.paramsD['regionid'])
                    
                    if layer.locus.locus != self.process.proc.paramsD['regionid']:
                        layer.locus.locus = self.process.proc.paramsD['regionid']
                    
                    if layer.locus.path != self.process.proc.paramsD['regionid']:
                        layer.locus.path = self.process.proc.paramsD['regionid']
                                                                      

                    layer._SetPath()
                    
                    if self.verbose:
                        print ('        filepath', layer.FPN)
   
                    layer.CreateVectorAttributeDef(fieldDD)
                    layer._SetBounds(query['epsg'],query['minlon'],query['minlat'],query['maxlon'],query['maxlat'])

                    projection = mj_gis.MjProj()
                    projection.SetFromEPSG(query['epsg'])

                    if not layer._Exists() or self.process.proc.overwrite:
                        mj_gis.CreateESRIPolygonPtL(layer.FPN, layer.fieldDefL, layer.BoundsPtL, projection.proj_cs, query['regionid'])          
                    boundsD = mj_gis.GetFeatureBounds(layer.FPN,'REGIONID')
                    
                    #Set lonlat projection
                    lonlatproj = mj_gis.MjProj()
                    lonlatproj.SetFromEPSG(4326)
                    
                    #Get the corners in lonlat
                    llD = mj_gis.ReprojectBounds(layer.BoundsPtL,projection.proj_cs,lonlatproj.proj_cs)

                    session._InsertDefRegion(self.process, layer, query, boundsD[query['regionid']], llD, self.process.proc.overwrite, self.process.proc.delete )

    def _DefaultRegFromVec(self,session):
        for locus in self.process.dstLayerD:
            if self.verbose:
                print ('locus',locus)
            for datum in self.process.dstLayerD[locus]:
                if self.verbose:
                    print ('datum',datum)
                for comp in self.process.dstLayerD[locus][datum]:

                    dstLayer = self.process.dstLayerD[locus][datum][comp]
                    #if not dstLayer._Exists() or self.process.proc.overwrite:
                    srcLayer = self.process.srcLayerD[locus][datum][comp]
                    if not path.isfile(srcLayer.FPN):
                        exitstr = 'setup_process_class: No source layer in _DefaultRegFromVec', srcLayer.FPN
                        exit(exitstr)
                    p = self.process.params
                    fieldL = [p.idcol, p.namecol, p.categorycol, p.parentidcol, p.parentcatcol, p.stratumcol,p.titlecol,p.labelcol]
                    fieldD = mj_gis.GetFeatureAttributeList(srcLayer.FPN, fieldL, p.idcol)
                    if not fieldD:
                        exit('setup_process_class: fieldD failed in _DefaultRegFromVec')
                    for key in fieldD:
                        fieldDD = self._SetfieldD( str(fieldD[key][p.idcol]), str(fieldD[key][p.namecol]), str(fieldD[key][p.categorycol]), int(fieldD[key][p.stratumcol]), str(fieldD[key][p.parentidcol]),str(fieldD[key][p.parentcatcol]) )
                        regionid, regioncat, parentid, regionname, parentcat = fieldD[key][p.idcol],fieldD[key][p.categorycol],fieldD[key][p.parentidcol],fieldD[key][p.namecol],fieldD[key][p.parentcatcol]
                        parentid = str(parentid.lower()).replace(' ', '-')
                        #regionname = str(regionname.lower()).replace(' ', '-')
                        regionid = str(regionid.lower()).replace(' ', '-')
                        if self.verbose:
                            print ('        forcing dst region from to',locus, regionid)

                        dstLayer.locus.locus = regionid
                        dstLayer.locus.path = regionid
                        if self.verbose:
                            print ('        forcing comp band and prefix region to "roi"')
                        #print ('dstLayer.comp.band',dstLayer.comp.band)
                        dstLayer.comp.band = dstLayer.comp.prefix = 'roi'
                        dstLayer.comp._SetCompid()
                        #print ('dstLayer.comp.compid',dstLayer.comp.compid)
                        dstLayer._SetPath()

                        dstLayer.CreateVectorAttributeDef(fieldDD)
                        fieldname = p.idcol
                        valueLL = [[fieldD[key][p.idcol]]]
                        if not dstLayer._Exists() or self.process.proc.overwrite:
                            mj_gis.ExtractFeaturesToNewDS(srcLayer.FPN, dstLayer.FPN,fieldname,valueLL, dstLayer.fieldDefL)
                            fieldname = 'REGIONID'
                            #Get the epsg and bounds
                            boundsD = mj_gis.GetFeatureBounds(dstLayer.FPN,fieldname) 
                            if len(boundsD) != 1:
                                exitstr = 'Default regions must consist on only one (1) feature (polygon or multipolygon): %s' %(dstLayer.FPN)
                                exit(exitstr)
                            projection = mj_gis.GetVectorProjection(dstLayer.FPN)
                            k = list(boundsD)[0]
                            bounds = boundsD[k]
    
                            dstLayer._SetBounds(projection.epsg,boundsD[k][0], boundsD[k][1], boundsD[k][2], boundsD[k][3] )  
                  
                            #Set lonlat projection
                            lonlatproj = mj_gis.MjProj()
                            lonlatproj.SetFromEPSG(4326)
                            #Get the corners in lonlat
                            llD = mj_gis.ReprojectBounds(dstLayer.BoundsPtL,projection.proj_cs,lonlatproj.proj_cs)
    
                            title = label = 'default region %s' %(regionid)
                            query = {'regionname':regionname,'regioncat':regioncat, 'parentid':parentid, 'parentcat':parentcat,'regionid':regionid, 'title':title,'label':label,'epsg':projection.epsg}
                            session._InsertDefRegion(self.process, dstLayer, query, bounds, llD,self.process.overwrite,self.process.delete )
                                                        
                        else:
                            boundsD = mj_gis.GetFeatureBounds(dstLayer.FPN,"REGIONID")
                            key = list(boundsD.keys())[0]
                            bounds = boundsD[key]
                            minx,miny,maxx,maxy = bounds
                            llparams = ['ullon','ullat','urlon','urlat','lrlon','lrlat','lllon','lllat']
                            llvalues = [minx,maxy,maxx,maxy,maxx,miny,minx,miny]
                            llD = dict(zip(llparams,llvalues))
                            title = label = 'default region %s' %(regionid)
                            query = {'regionname':regionname,'regioncat':regioncat, 'parentid':parentid, 'parentcat':parentcat,'regionid':regionid, 'title':title,'label':label,'epsg':4326}

                            session._InsertDefRegion(self.process, dstLayer, query, bounds, llD,self.process.overwrite,self.process.delete )

                            if self.verbose:
                                printstr = '    Layer %s already exists, skipping' %(dstLayer.FPN)
                                print (printstr)

    def _SetfieldD(self,regionid,regionname,regioncat,stratum,parentid,parentcat):
        #TGTODO SHOULD BE FROM DB
        fieldDD = {}
        fieldDD['REGIONID'] = {'name':'REGIONID', 'type':'string','width':32,'precision':0,'transfer':'constant','source':regionid }
        fieldDD['NAME'] = {'name':'NAME', 'type':'string','width':64,'precision':0,'transfer':'constant','source':regionname }
        fieldDD['CATEGORY'] = {'name':'CATEGORY', 'type':'string','width':32,'precision':0,'transfer':'constant','source':regioncat }
        fieldDD['STRATUM'] = {'name':'STRATUM', 'type':'integer','width':4,'precision':0,'transfer':'constant','source':stratum }
        fieldDD['PARENTID'] = {'name':'PARENTID', 'type':'string','width':32,'precision':0,'transfer':'constant','source':parentid }
        fieldDD['PARENTCAT'] = {'name':'PARENTCAT', 'type':'string','width':32,'precision':0,'transfer':'constant','source':parentcat }
        return fieldDD
          
          
class ProcessDefaultRegions:
    '''
    '''  
    def __init__(self, processparams ):
        
        self.processparams = processparams
        
        #direct to subprocess
        if self.processparams.processid == 'RegionCategories':
            
            self._InsertRegionCat(self.process)
            
        #elif self.process.proc.processid == 'defaultregion':
        elif self.processparams.processid == 'DefaultRegionFromCoords':
        
            self._DefaultRegionFromCoords()
            
        elif self.process.proc.processid == 'defaultregionfromvector':
            self._DefaultRegFromVec(session)
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
                    
    def _DefaultRegionFromCoords(self):
        '''
        '''
        
        locus = 'hej'
        for locus in self.process.dstLayerD:
            
            for datum in self.process.dstLayerD[locus]:
                
                for comp in self.process.dstLayerD[locus][datum]:
                    
                    query = self.process.proc.paramsD
                    fieldDD = self._SetfieldD( query['regionid'], query['regionname'], query['regioncat'], query['stratum'], query['parentid'], query['parentcat'])
                    layer = self.process.dstLayerD[locus][datum][comp]

                    #The destination region must be forced,this is because the locus to be created did not exists when chekcing for the feault locus
                    if self.verbose:
                        print ('        forcing region from to',locus, self.process.proc.paramsD['regionid'])
                    
                    if layer.locus.locus != self.process.proc.paramsD['regionid']:
                        layer.locus.locus = self.process.proc.paramsD['regionid']
                    
                    if layer.locus.path != self.process.proc.paramsD['regionid']:
                        layer.locus.path = self.process.proc.paramsD['regionid']
                                                                      

                    layer._SetPath()
                    
                    if self.verbose:
                        print ('        filepath', layer.FPN)
   
                    layer.CreateVectorAttributeDef(fieldDD)
                    layer._SetBounds(query['epsg'],query['minlon'],query['minlat'],query['maxlon'],query['maxlat'])

                    projection = mj_gis.MjProj()
                    projection.SetFromEPSG(query['epsg'])

                    if not layer._Exists() or self.process.proc.overwrite:
                        mj_gis.CreateESRIPolygonPtL(layer.FPN, layer.fieldDefL, layer.BoundsPtL, projection.proj_cs, query['regionid'])          
                    boundsD = mj_gis.GetFeatureBounds(layer.FPN,'REGIONID')
                    
                    #Set lonlat projection
                    lonlatproj = mj_gis.MjProj()
                    lonlatproj.SetFromEPSG(4326)
                    
                    #Get the corners in lonlat
                    llD = mj_gis.ReprojectBounds(layer.BoundsPtL,projection.proj_cs,lonlatproj.proj_cs)

                    session._InsertDefRegion(self.process, layer, query, boundsD[query['regionid']], llD, self.process.proc.overwrite, self.process.proc.delete )

    def _DefaultRegFromVec(self,session):
        for locus in self.process.dstLayerD:
            if self.verbose:
                print ('locus',locus)
            for datum in self.process.dstLayerD[locus]:
                if self.verbose:
                    print ('datum',datum)
                for comp in self.process.dstLayerD[locus][datum]:

                    dstLayer = self.process.dstLayerD[locus][datum][comp]
                    #if not dstLayer._Exists() or self.process.proc.overwrite:
                    srcLayer = self.process.srcLayerD[locus][datum][comp]
                    if not path.isfile(srcLayer.FPN):
                        exitstr = 'setup_process_class: No source layer in _DefaultRegFromVec', srcLayer.FPN
                        exit(exitstr)
                    p = self.process.params
                    fieldL = [p.idcol, p.namecol, p.categorycol, p.parentidcol, p.parentcatcol, p.stratumcol,p.titlecol,p.labelcol]
                    fieldD = mj_gis.GetFeatureAttributeList(srcLayer.FPN, fieldL, p.idcol)
                    if not fieldD:
                        exit('setup_process_class: fieldD failed in _DefaultRegFromVec')
                    for key in fieldD:
                        fieldDD = self._SetfieldD( str(fieldD[key][p.idcol]), str(fieldD[key][p.namecol]), str(fieldD[key][p.categorycol]), int(fieldD[key][p.stratumcol]), str(fieldD[key][p.parentidcol]),str(fieldD[key][p.parentcatcol]) )
                        regionid, regioncat, parentid, regionname, parentcat = fieldD[key][p.idcol],fieldD[key][p.categorycol],fieldD[key][p.parentidcol],fieldD[key][p.namecol],fieldD[key][p.parentcatcol]
                        parentid = str(parentid.lower()).replace(' ', '-')
                        #regionname = str(regionname.lower()).replace(' ', '-')
                        regionid = str(regionid.lower()).replace(' ', '-')
                        if self.verbose:
                            print ('        forcing dst region from to',locus, regionid)

                        dstLayer.locus.locus = regionid
                        dstLayer.locus.path = regionid
                        if self.verbose:
                            print ('        forcing comp band and prefix region to "roi"')
                        #print ('dstLayer.comp.band',dstLayer.comp.band)
                        dstLayer.comp.band = dstLayer.comp.prefix = 'roi'
                        dstLayer.comp._SetCompid()
                        #print ('dstLayer.comp.compid',dstLayer.comp.compid)
                        dstLayer._SetPath()

                        dstLayer.CreateVectorAttributeDef(fieldDD)
                        fieldname = p.idcol
                        valueLL = [[fieldD[key][p.idcol]]]
                        if not dstLayer._Exists() or self.process.proc.overwrite:
                            mj_gis.ExtractFeaturesToNewDS(srcLayer.FPN, dstLayer.FPN,fieldname,valueLL, dstLayer.fieldDefL)
                            fieldname = 'REGIONID'
                            #Get the epsg and bounds
                            boundsD = mj_gis.GetFeatureBounds(dstLayer.FPN,fieldname) 
                            if len(boundsD) != 1:
                                exitstr = 'Default regions must consist on only one (1) feature (polygon or multipolygon): %s' %(dstLayer.FPN)
                                exit(exitstr)
                            projection = mj_gis.GetVectorProjection(dstLayer.FPN)
                            k = list(boundsD)[0]
                            bounds = boundsD[k]
    
                            dstLayer._SetBounds(projection.epsg,boundsD[k][0], boundsD[k][1], boundsD[k][2], boundsD[k][3] )  
                  
                            #Set lonlat projection
                            lonlatproj = mj_gis.MjProj()
                            lonlatproj.SetFromEPSG(4326)
                            #Get the corners in lonlat
                            llD = mj_gis.ReprojectBounds(dstLayer.BoundsPtL,projection.proj_cs,lonlatproj.proj_cs)
    
                            title = label = 'default region %s' %(regionid)
                            query = {'regionname':regionname,'regioncat':regioncat, 'parentid':parentid, 'parentcat':parentcat,'regionid':regionid, 'title':title,'label':label,'epsg':projection.epsg}
                            session._InsertDefRegion(self.process, dstLayer, query, bounds, llD,self.process.overwrite,self.process.delete )
                                                        
                        else:
                            boundsD = mj_gis.GetFeatureBounds(dstLayer.FPN,"REGIONID")
                            key = list(boundsD.keys())[0]
                            bounds = boundsD[key]
                            minx,miny,maxx,maxy = bounds
                            llparams = ['ullon','ullat','urlon','urlat','lrlon','lrlat','lllon','lllat']
                            llvalues = [minx,maxy,maxx,maxy,maxx,miny,minx,miny]
                            llD = dict(zip(llparams,llvalues))
                            title = label = 'default region %s' %(regionid)
                            query = {'regionname':regionname,'regioncat':regioncat, 'parentid':parentid, 'parentcat':parentcat,'regionid':regionid, 'title':title,'label':label,'epsg':4326}

                            session._InsertDefRegion(self.process, dstLayer, query, bounds, llD,self.process.overwrite,self.process.delete )

                            if self.verbose:
                                printstr = '    Layer %s already exists, skipping' %(dstLayer.FPN)
                                print (printstr)

    def _SetfieldD(self,regionid,regionname,regioncat,stratum,parentid,parentcat):
        #TGTODO SHOULD BE FROM DB
        fieldDD = {}
        fieldDD['REGIONID'] = {'name':'REGIONID', 'type':'string','width':32,'precision':0,'transfer':'constant','source':regionid }
        fieldDD['NAME'] = {'name':'NAME', 'type':'string','width':64,'precision':0,'transfer':'constant','source':regionname }
        fieldDD['CATEGORY'] = {'name':'CATEGORY', 'type':'string','width':32,'precision':0,'transfer':'constant','source':regioncat }
        fieldDD['STRATUM'] = {'name':'STRATUM', 'type':'integer','width':4,'precision':0,'transfer':'constant','source':stratum }
        fieldDD['PARENTID'] = {'name':'PARENTID', 'type':'string','width':32,'precision':0,'transfer':'constant','source':parentid }
        fieldDD['PARENTCAT'] = {'name':'PARENTCAT', 'type':'string','width':32,'precision':0,'transfer':'constant','source':parentcat }
        return fieldDD
           
 
