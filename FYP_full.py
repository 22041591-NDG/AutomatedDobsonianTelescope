#!/usr/bin/python

from Adafruit_MotorHAT import Adafruit_MotorHAT, Adafruit_DCMotor, Adafruit_StepperMotor
#from MPU6050 import MPU6050
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u 
import time as tt
from astropy.time import Time
from datetime import *
from numpy import *
from math import *
pi = 3.1415926542

# coords in format [ RA, declination]
# RA format [hhmmss] h = hours
# dec format[ddmmss] d = degrees
userSavedCoord = ['******','000000','000000']
sirius     = ['064609','-164435']
alpha_cen  = ['144106','-605544']
rigel      = ['051538','-081021']
vega       = ['183742','384832']
betelgeuse = ['055624','072447']
canopus    = ['062428','-524205']
aldebaran  = ['043714','163320']
altair     = ['195153','085550']

def keywordsNames():
    print('The keywords availble on this system are:')
    print('\tSirius\n\tRigel\n\tVega\n\tBetelgeuse\n\tCanopus\n\tAldebaran\n\tAltair')

def readLocation():
    with open('latitude.txt','r') as file:
        lat = file.read()
        file.close()
    with open('longitude.txt','r') as file:
        long = file.read()
        file.close()
    return [float(lat), float(long)]

def writeLocation(lat,long):
    with open('latitude.txt','w') as file:
        file.write(str(lat))
        file.close()
    with open('longitude.txt','w') as file:
        file.write(str(long))
        file.close()

def getKeywordEC(ra_in, dec_in):
    ra  = hourToAngle(ra_in)
    dec = stringToAngle(dec_in)

def checkKeyword(key):
    if(key == userSavedCoord[0]):
        return[userSavedCoord[1], userSavedCoord[2]]
    elif(key.upper() == 'SIRIUS'):
        return sirius 
    elif(key.upper() == 'ALPHA_CENTAURI'):
        return alpha_cen
    elif(key.upper() == 'VEGA'):
        return vega
    elif(key.upper() == 'RIGEL'):
        return rigel
    elif(key.upper() == 'BETELGEUSE'):
        return betelgeuse
    elif(key.upper() == 'CANOPUS'):
        return canopus
    elif(key.upper() == 'ALDEBARAN'):
        return aldebaran
    elif(key.upper() == 'ALTAIR'):
        return altair

def stringToAngle(a):
    negPos = a[0]
    if(negPos == '-'):
        dd = a[1:3]
        mm = a[3:5]
        ss = a[5:7]
        angle = -(int(dd) +int(mm)/60+int(ss)/3600)
        return angle
    else:
        dd = a[0:2]
        mm = a[2:4]
        ss = a[4:6]
        angle = int(dd)+(int(mm)/60)+(int(ss)/3600)
        return angle

def hourToAngle(a):
    negPos = a[0]
    if(negPos == '-'):
        hh = a[1:3]
        mm = a[3:5]
        ss = a[5:7]
        angle = -(int(hh)*15 +int(mm)*0.25+int(ss)*0.004167)
        return angle
    else:
        hh = a[0:2]
        mm = a[2:4]
        ss = a[4:6]
        angle = int(hh)*15+(int(mm)*0.25)+(int(ss)*0.004167)
        return angle

def angleToOrient(az, alt):
    z = float(az /360)
    if(z>1):
        z = z % 1
    elif(z<-1):
        z = z % -1
    y = float(alt/90)
    if(y>1):
        y = y % 1
    if (y<-1):
        y = y%-1
    return[z, y]

lat  = -33.9328
latitude = -33.9328
longitude = 19

def ECtoHC(ra, dec):
    print('Start EC to HC')
    today = datetime.now()
    t     = Time(today, scale = 'utc', location = (lat ,longitude))
    a     = t.sidereal_time('apparent', longitude)
    LST   = a.degree - 30

    Hour = LST - ra
    eleRad = arcsin(sin(radians(lat))*sin(radians(dec)) + cos(radians(lat))*cos(radians(dec))*cos(radians(Hour)) )
    azRad  = arcsin( ( cos(radians(dec)) *sin(radians(Hour)) )/cos(eleRad)  )
    azRad2 = arctan2( sin(radians(Hour)) , sin(radians(lat))*cos(radians(Hour)) - cos(radians(lat))*tan(radians(dec)) )
    print('Done with calculations')
    if(azRad != azRad2):
        azRad = pi + azRad2
    elif(azRad == azRad2):
        azRad = 2*pi - azRad2
    alt = degrees(eleRad)
    az  = degrees(azRad)
    print('HC : az = ',az,'    el = ',alt)
    return [az, alt]

def ECtoHCt(ra, dec):
    sc = SkyCoord(ra,dec, frame = 'icrs', unit = 'deg')
    me_loc = EarthLocation(lat = latitude*u.deg, lon = longitude*u.deg, height = 93*u.m)
    otime  = Time.now()
    frame  = AltAz(obstime = otime, location = me_loc)
    ansHC  = sc.transform_to(frame)
    return [ansHC.az.deg, ansHC.alt.deg]

import atexit

# create a default object, no changes to I2C address or frequency
mh = Adafruit_MotorHAT()

# recommended for auto-disabling motors on shutdown!
def turnOffMotors():
    mh.getMotor(1).run(Adafruit_MotorHAT.RELEASE)
    mh.getMotor(2).run(Adafruit_MotorHAT.RELEASE)
    mh.getMotor(3).run(Adafruit_MotorHAT.RELEASE)
    mh.getMotor(4).run(Adafruit_MotorHAT.RELEASE)

atexit.register(turnOffMotors)

st1       = mh.getStepper(200, 1)
myStepper = mh.getStepper(200, 2)  # 200 steps/rev, motor port #1
myStepper.setSpeed(50)             # 40 RPM
st1.setSpeed(10)

def decimalCheck(a):
    len1 = len(a)
    pos = 0
    decimal = False
    while(1):
        if(a[pos-1] == '.'):
            decimal = True
        pos = pos + 1
        if(pos > len1):
            return decimal
            break

accel = [0]*3
gyro  = [0]*3
def getFeedback():
    cnt = 0
    maxVal = 0
    hor_orient = 0
    curr_time = time.time()
    while(True):
        prev_time = curr_time
        curr_time = time.time()
        time_diff = curr_time - prev_time
        accel = mpu.get_acceleration()
        gyro  = mpu.get_rotation()
        #print('x: ',gyro[0],'\ty: ',gyro[1],'\tz: ', gyro[2])
        hor_orient = hor_orient + gyro[2]*time_diff
        if(maxVal < accel[2]/16450):
            maxVal = accel[2]/16450   
        if(cnt > 2500):
            print('x: ',accel[0]/18156,'\ty: ',accel[1]/16450,'\tz: ', accel[2]/16450)
            cnt = 0
            print('max = ',maxVal)
            print('Orientation : ',hor_orient)
        else:
            print('max = ',maxVal)
            cnt = cnt+1

def numbersOnly(a):
    cnt = 0
    letterFlag = 0
    decimalFlag = 0
    while(cnt < len(a) ):
        numbers = ['1','2','3','4','5','6','7','8','9','0','.','-']
        num = 0
        prelimFlag = 0
        if(a[cnt] == '.'):
            decimalFlag = decimalFlag +1
        while(num < 12):
            if(a[cnt] != numbers[num]):
                prelimFlag = prelimFlag +1
            num = num +1
        if(prelimFlag == 12):
            letterFlag = letterFlag +1
        cnt = cnt +1
    if(letterFlag != 0):
        return False
    elif(decimalFlag > 1):
        return False
    else:
        return True

def getInput():    
    intype = raw_input('\nMENU\n1. Horizontal Coordinates\n2. Equatorial Coordinates\n3. Keywords\n4. Change Latitude and Longitude\n0. Exit Program\nMenu Number Selected: ')
    if( intype != '0' and intype != '1' and intype != '2' and intype != '3' and intype != '4'):
        return [0,0,'error1']
    if (int(intype) == 0 ):
        return[0,0,'exit']
    elif (int(intype) == 1 ):
        alt = (raw_input('Enter altitude: '))
        if(numbersOnly(alt) ):
	    alt = float(alt)
        else:
            print('A letter was entered.\nInvalid input.')
            return[0,0,'error2']
        az  = (raw_input('Enter azimuth: '))
        if(numbersOnly(az) ):
            az  = float(az)
        else:
            print('A letter was entered.\nInvalid input.')
            return[0,0,'error2']
	return[az , alt, 'HC']
    elif(int(intype) == 2):
        asc = raw_input('Enter ascension(hhmmss): ')
        dec = raw_input('Enter declination(ddmmss): ')
        if(len(dec)<6):
            print('dec too short')
            return[0,0,'error3']
        elif(len(asc)<6 ): 
            print('asc too short')
            return[0,0,'error3']
        elif(len(dec)>7 ): 
            print('dec too long')
            return[0,0,'error3']
        elif(len(asc)>7):
            print('asc too long')
            return[0,0,'error3']
        if(numbersOnly(asc) and numbersOnly(dec) ):
            return[asc, dec,'EC']
        else:
            print('A letter was entered.\nInvalid input.')
            return[0,0,'error2']
    elif(int(intype) == 3):
        keywordsNames()
        key = raw_input('Enter key word: ')
        return[1,key,'key']
    elif(int(intype) == 4):
        newLat  = raw_input('Input new latitude in degrees: ')
        newLong = raw_input('Input new longitude in degrees: ')
        number1 = numbersOnly(newLat)
        if(number1):
            decimal = decimalCheck(newLat)
            if(decimal):
                lat = float(newLat)
            else:
                lat = int(newLat)
        else:
            return [0,0,'error2']
        number2 = numbersOnly(newLong)
        if(number2):
            decimal = decimalCheck(newLong)
            if(decimal):
                longitude = float(newLong)
            else:
                longitude = int(newLong)
            print('new lat = ',lat,', new long = ',longitude)
        else:
            return [0,0,'error2']
        return [lat,longitude,'newLatLong']
    else:
        return [0,0,'error2']

hor = 0 #'''var for current horizontal angle'''
ver = 0 #'''var for current vertical angle'''

desired_ra  = 0 #'''var for current desired horizontal angle'''
desired_dec = 0 #'''var for current desired vertical angle'''

mode = 0        # 0 = HC ; 1 = EC ; 2 = keyword

def turnMotors(dy,dz,hor,ver):
    daz  = 4*(dy - hor)
    tt.sleep(2)
    dalt = dz - ver
    if(daz == 0 and dalt == 0):
        print('No movement')
    if(daz < 0):
        numSteps = math.ceil(-daz / 1.8)
#        print('turn HOR ',numSteps,' steps')
        hor = hor - (numSteps*1.8)/4
        myStepper.step(int(numSteps), Adafruit_MotorHAT.BACKWARD, Adafruit_MotorHAT.SINGLE)
    elif(daz > 0):
        numSteps = ceil(daz/1.8)
#        print('turn HOR ',numSteps,' steps')
        hor = hor + ( numSteps*1.8)/4
	myStepper.step(int(numSteps), Adafruit_MotorHAT.FORWARD, Adafruit_MotorHAT.SINGLE)
    tt.sleep(1)
    if(dalt < 0):
        numSteps = ceil(-dalt / 1.8)
        normDiff = ver - numSteps*1.8
        more     = numSteps + 1
        moreDiff = ver - more*1.8
        less     = numSteps - 1
        lessDiff = ver - less*1.8
        if( abs(moreDiff-dz)<abs(normDiff-dz) ):
            numSteps = more
        if( abs(lessDiff-dz)<abs(normDiff-dz) ):
            numSteps = less
#        print('turn VER ',numSteps,' steps')
        ver = ver - numSteps*1.8
	st1.step(int(numSteps),Adafruit_MotorHAT.FORWARD, Adafruit_MotorHAT.MICROSTEP)
    elif(dalt > 0):
        numSteps = ceil(dalt / 1.8)
        normDiff = ver + numSteps*1.8
        more     = numSteps + 1
        moreDiff = ver + more*1.8
        less     = numSteps - 1
        lessDiff = ver + less*1.8
        if( abs(moreDiff-dz)<abs(normDiff-dz) ):
            numSteps = more
        if( abs(lessDiff-dz)<abs(normDiff-dz) ):
            numSteps = less
#        print('turn VER ',numSteps,' steps')
        ver = ver + numSteps*1.8
	st1.step(int(numSteps), Adafruit_MotorHAT.BACKWARD, Adafruit_MotorHAT.MICROSTEP)
    print('New coords: hor = %.2f' % hor,'     ver = %.2f' % ver)
    turnOffMotors()
    return [hor, ver]

print('AUTOMATED DOBSONIAN MOUNT NEWTONIAN REFLECTOR TELESCOPE\n')
tt.sleep(1)
print('Please enter the menu number you want to input')

'''
print('Telescope Set to latitude of Stellenbosch University.')
print('Change is possible through setting.')
'''
time_then = 0

while(1):
    [lat,longitude] = readLocation()
    print('lat = ',lat,'  long = ',longitude)
    dz = 0
    dy = 0
    [hor1, ver1, inType] = getInput()
    if(inType == 'error1'):
        print('Inputs must be number of Menu Item')
    elif(inType == 'error2'):
        print('Error detected in input.Please try again.')
    elif(inType == 'error3'):
        print('Mistake in length of entry. Please try again.')
    elif(inType == 'newLatLong'):
        writeLocation(hor1, ver1)
    elif(inType == 'HC'):
        mode = 0
        az  = hor1
        alt = ver1
        dz = alt
        dy = az
        [hor, ver] = turnMotors(az, alt, hor, ver)
        time_then = tt.time()
    elif(inType == 'EC'):
        mode = 1
        ra = hourToAngle(hor1)
        dec = stringToAngle(ver1)
        [az, alt]   = ECtoHCt(ra, dec)
        [hor, ver]  = turnMotors(az, alt, hor, ver)
        print('Desired az = ',az,'    Desired alt = ',alt)
    elif(inType == 'key'):
        mode = 2
        [asc_1, dec_1]  = checkKeyword(ver1)
        ra  = hourToAngle(asc_1)
        dec = stringToAngle(dec_1)
        [az, alt]   = ECtoHCt(ra,dec)
        [hor, ver]  = turnMotors(az, alt, hor, ver)
        print('Desired az = ',az,'    Desired alt = ',alt)
    else:
        break
    time_now = tt.time()
    if((time_now - time_then > 60) and (mode != 0)):
         [az, alt]   = ECtoHCt(int(ra),int(dec))
#        [hor, ver]  = turnMotors(dy, dz, hor, ver)
    tt.sleep(0.1)

turnMotors(0,0,hor,ver)
print('Program Terminated')

'''
def canObjectBeSeen(ver):
    if(ver  <  0):
        print('Object cannot be viewed at this time of day.')
        print('Please try other coordinates.')
'''
