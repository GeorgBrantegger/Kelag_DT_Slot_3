                                     
# Testvolume                                            
# Depth of the whole structure is constant and given by the variable d
#                                           
#                                  
#        { x_1*d*h                                                      for       h <= h_1
# V(h) = { x_1*d*(h-h_1)+(x_2-x_1)*d*(h-h_1)**2/(2*(h_2-h_1) + V(h_1))  for h_1 < h <= h_2
#        { x_2*d*(h-h_2)+(x_3-x_2)*d*(h-h_2)**2/(2*(h_3-h_2) + V(h_2))  for h_2 < h <= h_3
#        { x_3*d*(h-h_3) + V(h_3)                                       for h_3 < h                           
#
#
#       { V/(x_1*d)                                                     for       V <= V_1
#h(V) = { (-b_2+sqrt(b_2**2-4*a_2*c_2)/(2*a_2))                         for V_1 < V <= V_2
#       { (-b_3+sqrt(b_3**2-4*a_3*c_3)/(2*a_3))                         for V_2 < V <= V_3
#       { (V-V_3)/(x_1*d)                                               for V_3 < V                           
#
# with
#     a_2 = 0.5*((x_2-x_1)*d)/(h_2-h_1)                                      
#     a_3 = 0.5*((x_3-x_2)*d)/(h_3-h_2)  
#                                           
#     b_2 = x_1*d-((x_2-x_1)*d*h_1)/(h_2-h_1)                                       
#     b_3 = x_2*d-((x_3-x_2)*d*h_2)/(h_3-h_2)                                       
#                                           
#     c_2 = ((x_2-x_1)*d*h_1**2)/(h_2-h_1)-h_1*x_1*d-(V-V_1)                                         
#     c_3 = ((x_3-x_2)*d*h_2**2)/(h_3-h_2)-h_2*x_2*d-(V-V_2)                                         
#                                           
#                                           
#                                           
#                                           
#                                           
#                    _____                      
#                   |     |                      |       
#                   |     |                      |       
#                   |     |                      | h_4 - h_3     
#                   |     |                     _|_      
#                 __| _ _ |__                    |      
#                /    x_3    \                   |     
#               /             \                  |     
#              /               \                 |    
#             /                 \                | h_3 - h_2  
#            /                   \               |    
#           /                     \              |     
#          /                       \             |     
#         /                         \            |      
#        /                           \          _|_      
#       <----------------------------->          | 
#        \            x_2            /           | h_2 - h_1  
#         \                         /            |  
#          \ _ _ _ _ _ _ _ _ _ _ _ /            _|_
#           |         x_1         |              | 
#           |                     |              | h_1  
#           |                     |              |   
#           |_____________________|             _|_





def test_1_parameters():
    h_1 = 10
    h_2 = 5     + h_1
    h_3 = 5     + h_2

    x_1 = 100
    x_2 = 101
    x_3 = 30

    d = 5

    vol_1 = x_1*d*h_1
    vol_2 = x_1*d*(h_2-h_1)+(x_2-x_1)*d*(h_2-h_1)**2/(2*(h_2-h_1)) + vol_1
    vol_3 = x_2*d*(h_3-h_2)+(x_3-x_2)*d*(h_3-h_2)**2/(2*(h_3-h_2)) + vol_2

    a_2 = 0.5*((x_2-x_1)*d)/(h_2-h_1)                                      
    a_3 = 0.5*((x_3-x_2)*d)/(h_3-h_2)  
                                            
    b_2 = x_1*d-((x_2-x_1)*d*h_1)/(h_2-h_1)                                       
    b_3 = x_2*d-((x_3-x_2)*d*h_2)/(h_3-h_2)                                       
                                            
    c_2 = ((x_2-x_1)*d*h_1**2)/(2*(h_2-h_1))-h_1*x_1*d                                        
    c_3 = ((x_3-x_2)*d*h_2**2)/(2*(h_3-h_2))-h_2*x_2*d

    return h_1,h_2,h_3,x_1,x_2,x_3,d,vol_1,vol_2,vol_3,a_2,a_3,b_2,b_3,c_2,c_3

def V_h_test_1(h):
    h_1,h_2,h_3,x_1,x_2,x_3,d,vol_1,vol_2,vol_3,a_2,a_3,b_2,b_3,c_2,c_3 = test_1_parameters()
    if h <= h_1:
        V = x_1*d*h
    elif (h_1 < h) and (h <= h_2):
        V = x_1*d*(h-h_1)+(x_2-x_1)*d*(h-h_1)**2/(2*(h_2-h_1)) + vol_1
    elif (h_2 < h) and (h <= h_3):
        V = x_2*d*(h-h_2)+(x_3-x_2)*d*(h-h_2)**2/(2*(h_3-h_2)) + vol_2
    elif (h_3 < h):
        V = x_3*d*(h-h_3) + vol_3

    return V

def h_V_test_1(V):
    h_1,h_2,h_3,x_1,x_2,x_3,d,vol_1,vol_2,vol_3,a_2,a_3,b_2,b_3,c_2,c_3 =test_1_parameters()
    if V <= vol_1:
        h = V/(x_1*d)
    elif (vol_1 < V) and (V <= vol_2):
        h = (-b_2+(b_2**2-4*a_2*(c_2-(V-vol_1)))**0.5)/(2*a_2) 
    elif (vol_2 < V) and (V <= vol_3):
        h = (-b_3+(b_3**2-4*a_3*(c_3-(V-vol_2)))**0.5)/(2*a_3)
    elif (vol_3 < V):
        h = (V-vol_3)/(x_3*d)+h_3
    return h


def test_2_parameters():
    x = 10
    d = 10
    return x,d

def V_h_test_2(h):
    x,d = test_2_parameters()
    return x*d*h

def h_V_test_2(V):
    x,d = test_2_parameters()
    return V/(x*d)

def show_parameters(test_version):
    h_1,h_2,h_3,x_1,x_2,x_3,d,vol_1,vol_2,vol_3,a_2,a_3,b_2,b_3,c_2,c_3 = test_1_parameters()
    x,d = test_2_parameters()

    if test_version == 1:
        print('h_1: ', h_1)
        print('h_2: ', h_2)
        print('h_3: ', h_3)
        print('x_1: ', x_1)
        print('x_2: ', x_2)
        print('x_3: ', x_3)
    elif test_version == 2:
        print('x: ', x)
        print('d: ', d)

    

