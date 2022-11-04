from MPU6050 import MPU6050
import time
from datetime import datetime

mpu = MPU6050()
mpu.dmp_initialize()
mpu.set_DMP_enabled(True)
mpu.int_status  = mpu.get_int_status()
packet_size     = mpu.DMP_get_FIFO_packet_size()
FIFO_count      = mpu.get_FIFO_count()
FIFO_buffer     = [0]*64
FIFO_count_list = list()


accel = [0]*3
gyro  = [0]*3
def setup():
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
        quat  = mpu.DMP_get_quaternion_int16(FIFO_buffer)
        grav  = mpu.DMP_get_gravity(quat)
        rpy   = mpu.DMP_get_euler_roll_pitch_yaw(quat, grav)
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
#            print('max = ',maxVal)
            cnt = cnt+1
        #time.sleep(2)

if __name__ == '__main__':
    print('Program is starting')
    setup()
    try:
        loop()
    except KeyboardInterrupt:
        pass
