from threading import Thread
from time import sleep
 
 
def myasync(f):
    def wrapper(*args, **kwargs):
        thr = Thread(target=f, args=args, kwargs=kwargs)
        thr.start()
    
    return wrapper
 
 
@myasync
def A():
    sleep(10)
    print("函数A睡了十秒钟。。。。。。")
    print("a function")
    
 
def B():
    print("b function")
    
    
A()
B()