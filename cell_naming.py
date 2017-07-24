import os

def quad_cell_lists():
    path = os.path.dirname(os.path.abspath(__file__))

    list_quads = []
    list_cryocells = [] 


    for point in range(1,9):
        for on in 'L R'.split():
            side = on + str(point)
            #~ print side
            with open(path+'/quad_cell_mapping/cell_naming_%s.txt'%side) as fid:
                lines = fid.readlines()
                
                for line in lines:
                    list_quads.append(line.split('\t')[0])
                    list_cryocells.append(line.split('\t')[1].replace('\n', ''))
                    
    return list_quads, list_cryocells 
                    

def cell_to_quad_dict():
    dic = {}
    list_quads, list_cryocells = quad_cell_lists()
    for cell, quad in zip(list_cryocells, list_quads):
        dic[cell] = quad
    return dic
    
def quad_to_cell_dict():
    dic = {}
    list_quads, list_cryocells = quad_cell_lists()
    for cell, quad in zip(list_cryocells, list_quads):
        dic[quad] = cell
    return dic
    

