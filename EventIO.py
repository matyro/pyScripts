
import struct
import binascii

import time

from enum import Enum

class Runh(Enum):
    header = 1
    run_number = 2
    date_of_run_begin = 3
    corsika_version = 4    
    num_of_obsv_level = 5
    obsv_level_height_i = 6 #0<=i<10 | add i
    slope_energy = 16
    energy_lower_limit = 17
    energy_upper_limit = 18
    flag_EGS4 = 19
    flag_NKG = 20
    cutoff_hadron = 21
    cutoff_muon = 22
    cutoff_electron = 23
    cutoff_photon = 24
    phys_constants_i = 25 #0<=i<50
    number_shower = 94
    xscatt = 248
    yscatt = 249 

class Evth(Enum):
    header = 1
    event_number = 2
    particle_id = 3
    total_energy = 4
    start_alt = 5
    number_of_first_target=6
    first_interaction_height=7
    px=8
    py=9
    pz=10
    zenith_angle=11
    azimuth_angle=12
    number_of_rndm_seq=13 #max 10
    seed_i=14             #0<= i < 10 |seed_5 = seed_i+3*i
    number_rndm_calls1=15    #mod 10^6
    number_rndm_calls2=16    #res    
    run_number=44
    date_of_run_begin=45
    corsika_version=46    
    number_of_obsv_level=47
    obsv_level_height_i=48 #0<=i<10 | obsv_level_height_2 = obsv_level_height_i + 2    
    slope_energy=58
    energy_lower_limit=59
    energy_upper_linit=60
    cutoff_hadron=61
    cutoff_muon=62
    cutoff_electron=63
    cutoff_photon=64

class Photon(Enum): #39 bunches 0<=n<=38 | pos + n*7
    num_cherenk_photon=1
    x_pos=2
    y_pos=3
    u_cos_x=4
    v_coy_y=5
    time_first_interaction=6
    prod_height=7




#sync-tag
#Item type (bits 0 to 15), a single bit (16) of user data, extension flag (bit 17), reserved bits (18 to 19), and version of this item type (bits 20 to 31).
#ident
#length
topLevelStr = struct.Struct('< i H h i i')
botLevelStr = struct.Struct('< H h i i')

CORHEADER_Str = struct.Struct('< 274f')
TELOFF_Str = struct.Struct('< I f 20f 20f')
PHOTON_Str = struct.Struct('< 8f') #put_real(tbunch.x, tbunch.y, tbunch.cx, tbunch.cy, tbunch.ctime, tbunch.zem, tbunch.photons, tbunch.lambda);
TELPOS_Str = struct.Struct('< I 4f')
TELOFF_Str = struct.Struct('< I f 20f 20f') #first all x, then all y 


headerDict = {"sync":0, "version":1, "version_extra":2, "id":3, "length":4}


idDict = {1200:"RUNH",
	  1201:"TELPOS",
          1202:"EVTH",
	  1203:"TELOFF",
          1204:"TELARRAY",
          1205:"PHOTONS",
          1206:"LAYOUT",
          1207:"TRIGTIME",
          1208:"PE",
          1209:"EVTE",
          1210:"RUNE",
          1211:"LONGI",
          1212:"INPUTCFG",
          1213:"TELARRAY_HEAD",
          1214:"TELARRAY_END",
          1215:"EXTRA_PARAM"}

idDict_inv = {v: k for k, v in idDict.items()}




def printHeader(head, sync=True):
    if sync:
        print("Sync: ", head[0], "type: ", head[1], "version: ", head[2], "ID: ", head[3], "length: ", head[4] & 0x3FFFFFFF)
        print("Block: ", idDict[head[1]])
        print("Only Subitems: ", (head[4] & 0b01000000000000000000000000000000)!=0 )	
        print("Ext. Field active: ",(head[2] & 0x20000) != 0)
    else:
        print("\tType: ", head[0], "Version:", head[1], "ID: ", head[2], "length: ", head[3] & 0x3FFFFFFF)
        print("\tBlock: ", idDict[head[0]])
        print("\tOnly Subitems: ", (head[3] & 0b01000000000000000000000000000000)!=0 )
        print("\tExt. Field active: ",(head[1] & 0x2000) != 0)



def readSubBlock(dataBytes, maxLength):
    subHeader = botLevelStr.unpack( dataBytes[0:12] )
    #printHeader(subHeader, False)
   
    data = dataBytes[12:(subHeader[3] & 0x3FFFFFFF)+12]

    subHeader = (0,)+subHeader
    return {"header":subHeader, "data":data}



def readBlock(filePtr):
    headerByte = filePtr.read( 16 )
    if len(headerByte) == 0:
        return {"header":(0,0,0,0), "data":b''}

    header = topLevelStr.unpack( headerByte )
    #printHeader(header, True)

    data = filePtr.read(header[4] & 0x3FFFFFFF)
    return {"header":header, "data":data}
    




def main():
    with open("cer000000", "rb") as f:            
        
        start = time.time()
        telarray_counter = 0
        event_counter = 0
 
        runh_data = ()
        evth_data = ()

        while 1:
            #print("\n")
            
            block = readBlock(f)

            if block["data"] == b'':
                break

            
            if (block["header"][4]&0b01000000000000000000000000000000)!=0:
                subBlock = readSubBlock(block["data"], block["header"][4] & 0x3FFFFFFF)
                if subBlock["header"][1] == idDict_inv["PHOTONS"]:                    
                    photonData = PHOTON_Str.iter_unpack(subBlock["data"][12:])                    
                   
            

            if block["header"][1] == idDict_inv["RUNH"]:
                runh_data = CORHEADER_Str.unpack(block["data"])
                #print(runh_data)

            if block["header"][1] == idDict_inv["EVTH"]:
                evth_data = CORHEADER_Str.unpack(block["data"])
                #print(evth_data)
                telarray_counter = 0
                event_counter = event_counter + 1

            if block["header"][1] == idDict_inv["TELOFF"]:
                teloff_data = TELOFF_Str.unpack(block["data"])
                #print(teloff_data)

            if block["header"][1] == idDict_inv["TELARRAY"]:
                telarray_counter = telarray_counter + 1
                #print(telarray_counter)
            
            if block["header"][1] == idDict_inv["TELPOS"]:
                telposs_data = TELPOS_Str.unpack(block["data"])
                #print(telposs_data)

           
            if event_counter % 100 == 0 and event_counter > 0:
                print("Average Events/s: ", event_counter/(time.time() - start))


        print("All Top-Level Events processed: ", event_counter)
        print("Reuse: ", evth_data[98])  
 
       

if __name__ == "__main__":
    main()
