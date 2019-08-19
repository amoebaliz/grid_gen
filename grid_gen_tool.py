from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import inv
from math import sin, cos, sqrt, atan2, asin, acos, atan, degrees, radians 
import time
class Click():
    def __init__(self, ax, func, button=1):
        self.ax=ax
        self.func=func
        self.button=button
        self.press=False
        self.move = False
        self.c1=self.ax.figure.canvas.mpl_connect('button_press_event', self.onpress)
        self.c2=self.ax.figure.canvas.mpl_connect('button_release_event', self.onrelease)
        self.c3=self.ax.figure.canvas.mpl_connect('motion_notify_event', self.onmove)

    def onclick(self,event):
        if event.inaxes == self.ax:
            if event.button == self.button:
                self.func(event)
    def onpress(self,event):
        self.press=True
    def onmove(self,event):
        if self.press:
            self.move=True
    def onrelease(self,event):
        if self.press and not self.move:
            self.onclick(event)
        self.press=False; self.move=False

def calc_dist(glon_A1,glat_A1,glon_B1,glat_B1):
    # Distance Calculation
    rlat_A1 = radians(glat_A1)
    rlon_A1 = radians(glon_A1)
    rlat_B1 = radians(glat_B1)
    rlon_B1 = radians(glon_B1)

    d_rlon = rlon_B1 - rlon_A1
    d_rlat = rlat_B1 - rlat_A1

    # Haversine formulation
    a = sin(d_rlat / 2)**2 + cos(rlat_A1) * cos(rlat_B1) * sin(d_rlon / 2)**2
    # c = 2 * asin2(sqrt(a)) = same as below based on pythagorean relationships
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    distance = Re * c

    # Calculate rotation matrices

    # Rotates second point to Prime Meridian
    R1 = np.array((( cos(rlon_A1), sin(rlon_A1), 0),\
                   (-sin(rlon_A1), cos(rlon_A1), 0),\
                   ( 0           , 0           , 1)))

    # Rotates second point to equator
    R2 = np.array(((cos(rlat_A1),  0, sin(rlat_A1)),\
                   (0           ,  1, 0           ),\
                   (sin(-rlat_A1), 0, cos(rlat_A1))))

    # Original 3D Coordinates of points
    A1 = coord_3D(rlon_A1, rlat_A1)
        
    B1 = coord_3D(rlon_B1, rlat_B1)
    # once-rotated 3D Coordinates of points
    A2 = np.dot(R1,A1) 
    
    B2 = np.dot(R1,B1) 
    # Lat/Lon of twice-rotated first point       
    rlat_A2, rlon_A2 = new_latlon(A2[0],A2[1],A2[2])

    rlat_B2, rlon_B2 = new_latlon(B2[0],B2[1],B2[2])

    #m.plot(degrees(rlon_A2),degrees(rlat_A2),'mo')
    #m.plot(degrees(rlon_B2),degrees(rlat_B2),'ro')
    #plt.draw()      
    # Twice-rotated 3D Coordinates of first point
    A2 = np.dot(R2,np.dot(R1,A1)) 
       
    B2 = np.dot(R2,np.dot(R1,B1)) 
    # Lat/Lon of twice-rotated first point       
    rlat_A2, rlon_A2 = new_latlon(A2[0],A2[1],A2[2])

    rlat_B2, rlon_B2 = new_latlon(B2[0],B2[1],B2[2])

    #m.plot(degrees(rlon_A2),degrees(rlat_A2),'mo')
    #m.plot(degrees(rlon_B2),degrees(rlat_B2),'ro')
    #plt.draw()
    #plt.show()

    # Calculate angle theta between A2 and equator
    d = (cos(rlat_B2)-cos(rlon_B2)*cos(distance/Re))/ \
        (sin(rlon_B2)*sin(distance/Re))
    if glat_B1-glat_A1 == 0:
       x_sign=1
    else: 
       x_sign=np.sign(glat_B1-glat_A1)

    theta = x_sign*acos(d)
    print theta   
    # Rotates first point to equator
    R3 = np.array(((1,              0,              0),\
                   (0,  np.cos(theta),  np.sin(theta)),\
                   (0, np.sin(-theta), np.cos(theta)))) 

    #A3 = np.dot(R3,np.dot(R2,np.dot(R1,A1)))

    #B3 = np.dot(R3,np.dot(R2,np.dot(R1,B1)))
    #rlat_A3, rlon_A3 = new_latlon(A3[0],A3[1],A3[2])

    #rlat_B3, rlon_B3 = new_latlon(B3[0],B3[1],B3[2])

    #m.plot(degrees(rlon_A3),degrees(rlat_A3),'bo',latlon=True)
    #m.plot(degrees(rlon_B3),degrees(rlat_B3),'bo',latlon=True)
    #plt.draw()

    # BACK ROTATION CHECK
    #B2 = np.dot(inv(R3),B3)

    #A2 = np.dot(inv(R3),A3) 
    #rlat_A2, rlon_A2 = new_latlon(A2[0],A2[1],A2[2])

    #rlat_B2, rlon_B2 = new_latlon(B2[0],B2[1],B2[2])

    #m.plot(degrees(rlon_A2),degrees(rlat_A2),'mo',latlon=True)
    #m.plot(degrees(rlon_B2),degrees(rlat_B2),'mo',latlon=True)
    #plt.draw()

    #B1 = np.dot(inv(R2),B2) 

    #A1 = np.dot(inv(R2),A2) 
    #rlat_A1, rlon_A1 = new_latlon(A1[0],A1[1],A1[2])

    #rlat_B1, rlon_B1 = new_latlon(B1[0],B1[1],B1[2])

    #m.plot(degrees(rlon_A1),degrees(rlat_A1),'ko',latlon=True)
    #m.plot(degrees(rlon_B1),degrees(rlat_B1),'ko',latlon=True)
    #plt.draw()

    #B0 = np.dot(inv(R1),B1) 

    #A0 = np.dot(inv(R1),A1) 
    #rlat_A0, rlon_A0 = new_latlon(A0[0],A0[1],A0[2])

    #rlat_B0, rlon_B0 = new_latlon(B0[0],B0[1],B0[2])

    #m.plot(degrees(rlon_A0),degrees(rlat_A0),'go',latlon=True)
    #m.plot(degrees(rlon_B0),degrees(rlat_B0),'go',latlon=True)
    #plt.draw()

        
    return distance, R1, R2, R3 

def grid_gen(event):
    global MEEP, glat_A1, glon_A1, rlons_3, R1, R2, R3
    c_glon, c_glat = m(event.xdata,event.ydata,inverse=True)
    # PLOT FIRST POINT A
    if MEEP<1:
       glon_A1 = c_glon; glat_A1 = c_glat
       #glon_A1 = event.xdata; glat_A1 = event.ydata
       m.plot(glon_A1,glat_A1,'bo',latlon=True)
       plt.draw()
       MEEP+=1

    # CALCULATE AND PLOT SECOND POINT (B) ALONG GC LINE
    #   1) Calculate distance of AB 
    #	2) Generate rotation matrices (R1, R2) to put AB 
    #      on equator with terminus at the prime meridian 
    #	   (becomes A'B')
    #	3) Calculate grid positions along A'B' given resolution and 
    #      starting at (0,0); ceil(dist(AB)/resolution). 
    #   4) Plot R1-1 R2-1 A'B' and new points
  
    elif MEEP<2:
       # Generate Rotation Matrices Based on lat/lons
       
       #m.drawgreatcircle(glon_A1, glat_A1, event.xdata, event.ydata,del_s=10,color='k', lw=2.)
       init_distance, R1, R2, R3 = calc_dist(glon_A1,glat_A1,c_glon,c_glat)

       npts = int(abs(np.ceil(init_distance/grid_res)))+1 

       # Longitudes (radians) are regularly spaced (grid_res) going
       # away from the Prime Meridian in the direction of x_sign
       rlons_3 = (grid_res/Re)*np.arange(npts)

       # Along equator all points are at lat = 0

       # Get x,y,z coordinates
       X3 = coord_3D(rlons_3, np.zeros(npts))
       # Rotate back to original position
       X1 = np.dot(inv(R3),X3)
       X1 = np.dot(inv(R2),np.dot(inv(R3),X3))
       X1 = np.dot(inv(R1),np.dot(inv(R2),np.dot(inv(R3),X3)))
       # Convert to geographic coordinates and plot
       glats = np.zeros(npts) 
       glons = np.zeros(npts)
       for nt in range(npts):
           lat, lon = new_latlon(X1[0,nt], X1[1,nt], X1[2,nt])

           glats[nt] = degrees(lat)
           glons[nt] = degrees(lon)

           map_x, map_y = m(degrees(lon),degrees(lat))
           if m.is_land(map_x,map_y):
              m.plot(glons[nt],glats[nt],'go',zorder=3,latlon=True)
           else:
              m.plot(glons[nt],glats[nt],'bo',zorder=3,latlon=True)

       m.drawgreatcircle(glons[0], glats[0], glons[-1], glats[-1],del_s=10,color='k', lw=2.)
       plt.draw()

       # Rotate back to original position
       X1 = np.dot(inv(R1),np.dot(inv(R2),np.dot(inv(R3),X3)))

       # Convert to geographic coordinates and plot
       glats = np.zeros(npts) 
       glons = np.zeros(npts)
       for nt in range(npts):
           lat, lon = new_latlon(X1[0,nt], X1[1,nt], X1[2,nt])
 
           glats[nt] = degrees(lat)
           glons[nt] = degrees(lon)

           if m.is_land(glons[nt],glats[nt]):
              m.plot(glons[nt],glats[nt],'go',zorder=3)
           else:
              m.plot(glons[nt],glats[nt],'bo',zorder=3)

       m.drawgreatcircle(glons[0], glats[0], glons[-1], glats[-1],del_s=10,color='k', lw=2.)
       plt.draw()
       MEEP+=1

    elif MEEP<3:
         C1 = coord_3D(radians(c_glon),radians(c_glat))
         C3 = np.dot(R3,np.dot(R2,np.dot(R1,C1)))

         # Determine direction to send meridians based on position 
         # of C relative to AB
         merid_dir = np.sign(C3[-1])

         rlat_C3, rlon_C3 = new_latlon(C3[0],C3[1],C3[2])

         # minimmum distance from equator:
         min_dist = Re*rlat_C3
         npts = int(abs(np.ceil(min_dist/grid_res)))+1
 
         rlons_3 = np.tile(rlons_3,(npts,1))
         rlats_3 = merid_dir*np.tile((grid_res/Re)*np.arange(npts).reshape(npts,1),(1,rlons_3.shape[1]))
         X3 = coord_3D(rlons_3, rlats_3)

         nj = X3.shape[1]
         ni = X3.shape[2]

         X3 = X3.reshape(3,nj*ni)
         # Rotate back to original position
         X1 = np.dot(inv(R1),np.dot(inv(R2),np.dot(inv(R3),X3)))
         X1 = X1.reshape(3,nj,ni)

         glats = np.zeros((nj,ni))
         glons = np.zeros((nj,ni))

         for jj in range(nj):
             for ii in range(ni):
                 lat, lon = new_latlon(X1[0,jj,ii], X1[1,jj,ii], X1[2,jj,ii])

                 glats[jj,ii] = degrees(lat)
                 glons[jj,ii] = degrees(lon)

         m.plot(glons,glats,'bo',latlon=True)
         plt.draw()
         plt.show()
         # SAVE AS 2D ARRAY... pre land-masking  
         #np.savez('test_grid', glats, glons)
         MEEP+=1

def coord_3D(rlon,rlat):
    
    x = Re*np.cos(rlat)*np.cos(rlon) 
    y = Re*np.cos(rlat)*np.sin(rlon)
    z = Re*np.sin(rlat) 
   
    return np.array(((x), (y), (z)))

def new_latlon(x,y,z):
    lat = asin(z/np.abs(Re))
    lon = atan(y/x)

    return lat,lon

MEEP=0
# approximate radius of earth in km
Re = 6378.13
 
fig = plt.figure(figsize=(11.7,8.3))
plt.subplots_adjust(left=0.05,right=0.95,top=0.90,bottom=0.05,wspace=0.15,hspace=0.05)
ax = plt.subplot(111)
m = Basemap(projection='npstere',boundinglat=-30,lon_0=270,resolution='l')
m.drawmapboundary(fill_color='azure')
m.fillcontinents(color='palegoldenrod',lake_color='azure')
m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,30.),labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,420.,60.),labels=[0,0,0,1])
grid_res = input("Grid resolution in km: ") 

click = Click(ax, grid_gen)
#fig.canvas.mpl_connect("button_press_event", onclick)
plt.show()
