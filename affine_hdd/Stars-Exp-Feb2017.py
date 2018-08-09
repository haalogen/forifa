
# coding: utf-8

# In[12]:

import os
import sys
from PIL import Image, ImageDraw
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import norm, pinv
# get_ipython().magic('matplotlib inline')
np.set_printoptions(precision=12, suppress=True)

# Command Line Input

CENTER_RAD = int(sys.argv[1]) # radius(px) of central part
print('CENTER_RAD:', CENTER_RAD)

# In[13]:

zoomed_coords = True # whether to divide coordinates by SCALE_FAC or not (zoomed coords or not)
center_only = True # use only central stars
SCALE_FAC = 4.0 # Scale factor of coordinates


# In[14]:

"""
Load star coords from txt-files
""";


# In[15]:

coords_folder = 'data/star_coords/2016nov-11_txt/' # folder with coords files
images_folder = 'data/stars/2016nov-11/'
results_folder = 'results/2016nov-11/exp-Feb-2017'


# In[16]:

# 2016nov-11 jpg
fnames = [
    "20161122-191517-359.txt",
    "20161122-201517-375.txt",
    "20161122-211517-375.txt",
    "20161122-221517-375.txt"
]


# In[17]:

date = fnames[0][:-4]


# In[18]:

im = Image.open(images_folder + "mod_" + date + "-1.jpg")
w, h = im.size
print("Image size:", w, h)




# In[20]:

x_c, y_c = w/2, h/2 # central point of img
print('x_c, y_c:', x_c, y_c)


# "1." Берем область (круг) с центром в середине кадра, разных радиусов.
# Считаем количество звезд. Определяем координаты выделенных звезд,
# стараясь брать звезды малого размера, занимающих один пиксель.
# Число звезд в кругах должно  отличаться на 3-4.
# 
# * R1 = 576px; 11 звезд учитыв-ся при подсчете коэфф
# * R2 = 476px; 8 звезд учитыв-ся при подсчете коэфф

# In[21]:

"""
Choose central stars (in central circle)
""";


# In[23]:

coords_list = []
for fname in fnames:
    piece = np.loadtxt(coords_folder + os.sep + fname)
    coords_list.append(piece)

coords = np.vstack(coords_list)
print('Zoomed In Star coordinates pairs (first 5):\n', coords[:5], '\n')

if zoomed_coords:
    coords /= float(SCALE_FAC)
    coords = coords.round()
    print('Normal Star coordinates pairs (first 5):\n', coords[:5], '\n')


if center_only:
    coords_center = []
    
    for i in range(coords.shape[0]):
        _lx = coords[i, 0]
        _ly = coords[i, 1]
        _rx = coords[i, 2]
        _ry = coords[i, 3]
        if         (_lx - x_c)**2 + (_ly - y_c)**2 <= CENTER_RAD**2 and         (_rx - x_c)**2 + (_ry - y_c)**2 <= CENTER_RAD**2:
            coords_center.append(coords[i])
    
    coords = np.vstack(coords_center)
    print('Normal Star coordinates pairs in center:\n', coords, '\n')
    print('Num of star pairs chosen:', len(coords))


# In[24]:

lX = coords[:, 0] # leftX coordinates
lY = coords[:, 1] # leftY coordinates
rX = coords[:, 2] # rightX coordinates
rY = coords[:, 3] # rightY coordinates

N = coords.shape[0] # number of pairs of points
M = coords.shape[1] # lX, lY, rX, rY == 4
print('Number of Star coordinates pairs:', N)


# In[25]:

# Plot chosen star pairs
scatter_original = Image.new('RGB', (w, h), 'lightgray')
s_pix = scatter_original.load()

ELL_RAD = 3

scx_c = scatter_original.width / 2
scy_c = scatter_original.height / 2
draw = ImageDraw.Draw(scatter_original)
draw.ellipse((scx_c-ELL_RAD, scy_c-ELL_RAD, scx_c+ELL_RAD, scy_c+ELL_RAD), fill='darkgreen')


# In[27]:

for i in range(N): # draw star points
    draw.ellipse((lX[i]-ELL_RAD, lY[i]-ELL_RAD, lX[i]+ELL_RAD, lY[i]+ELL_RAD), fill='red')
    draw.ellipse((rX[i]-ELL_RAD, rY[i]-ELL_RAD, rX[i]+ELL_RAD, rY[i]+ELL_RAD), fill='blue')

# Draw central part boundary
draw.ellipse((scx_c-CENTER_RAD, scy_c-CENTER_RAD, scx_c+CENTER_RAD, scy_c+CENTER_RAD), outline='black')


# In[28]:

# Original stars' position
# scatter_original


# "2." Для каждого из набора звезд в выделенных кругах вычисляем
# коэффициенты аффинного преобразования, коэффициенты
# преобразования аффинного+дисторсии 3 порядка и коэффициенты
# преобразования аффинного+дисторсии 3 порядка+дисторсии 5
# порядка. Формулу для дисторсии 5 порядка можно посмотреть,
# например, на сайте
# https://ru.wikipedia.org/wiki/%D0%94%D0%B8%D1%81%D1%82%D0%BE%D1%80%D1%81%D0%B8%D1%8F

# Для дисторсии высших порядков  ($F_{n} при \; n>3$) в формулу добавляют по одному слагаемому на каждую дисторсию нечётного порядка ($F_{3}, F_{5}, F_{7} \;$ и т.п.):
# 
# $
# \vec{R} = b_{0}\vec {r} + F_{3} r^{2} \vec{r} + F_{5} r^{4} \vec{r} + F_{7} r^{6} \vec{r} + \dots
# $

# и написать по аналогии с третьим порядком. Число звезд должно быть
# больше, чем половина числа неизвестных параметров.

# Вычисление проводится так.
# Запишем схему измерения
# $$
# \xi = Af + \nu ,
# $$

# задаем матрицу А, матрицу $\Sigma$, являющуюся диагональной, и на ее
# диагонали стоят дисперсии координат вектора погрешности $\nu$ . Здесь $f$ -- вектор неизвестных коэффициентов преобразований. Считаем эти дисперсии равными единице, т.е. $\Sigma = I$ -- единичная матрица

# Оценка коэффициентов вычисляется по формуле
# $$
#     \hat{f} = A^- \xi
# $$

# In[29]:

xi = np.zeros(2*N)

for i in range(N): # fill the xi vector
    xi[2*i] = rX[i]
    xi[2*i + 1] = rY[i]

print('xi:\n', xi)


# In[30]:

"""
Calculate coeff-s

a) Affine (at least 3 stars)
""";


# In[31]:

k = 6 # num of coeff-s

z = np.zeros(k)
arr = np.zeros((2*N, k)) # matrix A


# In[32]:

for i in range(N): # fill the A matrix
    
    arr[2*i] = [lX[i], lY[i], 0, 0, 1, 0]

    arr[2*i + 1] = [0, 0, lX[i], lY[i], 0, 1]


# In[33]:

np.set_printoptions(precision=2, suppress=True)
print('A:\n', arr[:4], '\n')
np.set_printoptions(precision=12, suppress=True)


# In[34]:

p_arr = pinv(arr)
z = np.dot(p_arr, xi)
print("""
Affine coefficients:
%.4f %.4f %.4f %.4f 
%.2f %.2f""" % tuple(z))
print('cond(A): ', np.linalg.cond(arr))


# In[ ]:




# In[ ]:




# In[35]:

"""
b) Affine + Ditortion 3rd order 
  (at least 4 stars)
""";


# In[36]:

k3 = 8

z3 = np.zeros(k3)
arr3 = np.zeros((2*N, k3)) # matrix A


# In[37]:

for i in range(N): # fill the A matrix
    dist_l = (lX[i]-x_c)**2 + (lY[i]-y_c)**2
    dist_r = (rX[i]-x_c)**2 + (rY[i]-y_c)**2

    zx1 = (lX[i] - x_c) * dist_l
    zx2 = (rX[i] - x_c) * dist_r
    arr3[2*i] = [lX[i], lY[i], 0, 0, 1, 0, -zx1, zx2]

    zy1 = (lY[i] - y_c) * dist_l
    zy2 = (rY[i] - y_c) * dist_r
    arr3[2*i + 1] = [0, 0, lX[i], lY[i], 0, 1, -zy1, zy2]


# In[39]:

np.set_printoptions(precision=2)
print('A:\n', arr3[:4], '\n')
np.set_printoptions(precision=4)


# In[41]:

p_arr3 = pinv(arr3)
z3 = np.dot(p_arr3, xi)
print("""
Affine coefficients + Ditortion 3rd order: 
%.4f %.4f %.4f %.4f 
%.2f %.2f 
%.2e %.2e""" % tuple(z3))
print('cond(A): ', np.linalg.cond(arr3))


# In[ ]:




# In[ ]:




# In[42]:

"""
c) Affine + Ditortion 3rd, 5th orders 
  (at least 5 stars)
""";


# In[43]:

k35 = 10

z35 = np.zeros(k35)
arr35 = np.zeros((2*N, k35)) # matrix A


# In[44]:

for i in range(N): # fill the A matrix
    dist_l = (lX[i]-x_c)**2 + (lY[i]-y_c)**2
    dist_r = (rX[i]-x_c)**2 + (rY[i]-y_c)**2

    zx1 = (lX[i] - x_c) * dist_l
    zx2 = (rX[i] - x_c) * dist_r
    wx1 = (lX[i] - x_c) * dist_l**2
    wx2 = (rX[i] - x_c) * dist_r**2

    arr35[2*i] = [lX[i], lY[i], 0, 0, 1, 0, -zx1, zx2, -wx1, wx2]

    zy1 = (lY[i] - y_c) * dist_l
    zy2 = (rY[i] - y_c) * dist_r
    wy1 = (lY[i] - y_c) * dist_l**2
    wy2 = (rY[i] - y_c) * dist_r**2

    arr35[2*i + 1] = [0, 0, lX[i], lY[i], 0, 1, -zy1, zy2, -wy1, wy2]


# In[45]:

np.set_printoptions(precision=2, suppress=True)
print('A:\n', arr35[:4], '\n')


# In[46]:

p_arr35 = pinv(arr35, rcond=1e-18)
z35 = np.dot(p_arr35, xi)
np.set_printoptions(precision=4)
print("""
Affine coefficients + Ditortion 3rd, 5th orders:

%.4f %.4f %.4f %.4f 
%.2f %.2f 
%.2e %.2e 
%.2e %.2e""" % tuple(z35))
print('cond(A): ', np.linalg.cond(arr35))


# In[ ]:




# In[ ]:




# Погрешность оценки коэффициентов определяется матрицей
# $$
#     \Sigma_f = \left( A^* \Sigma A \right)^{-1} 
#              = \left(A^* A \right)^{-1}
# $$

# In[47]:

"""
Calculate error matrix \Sigma_f
""";


# In[48]:

sig_mat = np.linalg.inv(arr.T.dot(arr))
# plt.matshow(sig_mat, cmap=plt.cm.gray)
np.set_printoptions(suppress=False)
print(sig_mat[:4, :4])


# In[ ]:




# In[49]:

sig_mat3 = np.linalg.inv(arr3.T.dot(arr3))
# plt.matshow(sig_mat3, cmap=plt.cm.gray)
np.set_printoptions(suppress=False)
print(sig_mat3[:4, :4])


# In[ ]:




# In[50]:

sig_mat35 = np.linalg.inv(arr35.T.dot(arr35))
# plt.matshow(sig_mat35, cmap=plt.cm.gray)
np.set_printoptions(suppress=False)
print(sig_mat35[:4, :4])


# In[ ]:




# Надежность модели определяется значением случайной величины
# $$
#     \tau = \| \xi - AA^{-} \xi \|^2
# $$

# In[51]:

"""
Calculate reliability of model \tau,
(2N - k) / tau,
N -- number of central stars 
k -- number of coefficients(6, 8 or 10)
""";


# In[52]:

N


# In[53]:

tau = np.linalg.norm(xi - arr.dot(p_arr).dot(xi) ) ** 2
print("tau:", tau)
print("(2N - k) / tau:", (2*N - k) / tau)


# In[ ]:




# In[54]:

tau3 = np.linalg.norm(xi - arr3.dot(p_arr3).dot(xi) ) ** 2
print("tau3:", tau3)
print("(2N - k3) / tau3:", (2*N - k3) / tau3)


# In[ ]:




# In[55]:

tau35 = np.linalg.norm(xi - arr35.dot(p_arr35).dot(xi) ) ** 2
print("tau35:", tau35)
print("(2N - k35) / tau35:", (2*N - k35) / tau35)


# In[ ]:




# Выходные данные программы:
# 
# 1) N – число звезд в круге
# 
# 2) оценка коэффициентов $\; \hat{f} = A^- \xi$
# 
# 3) матрица ошибок $\; \Sigma_f$
# 
# 4) значение $\; \tau$ , а лучше – значение $\; \frac{2N - k}{\tau}$,
# где k -- число коэффициентов (координат вектора f).

# In[ ]:




# "3." Вычисляем преобразованное изображение одного из кадров и
# накладываем его на изображение второго кадра. Для каждой звезды
# вычисляем разность координат (по x и по y) двух (совмещенных) изображений звезд.

# In[56]:

def affine_transform_point(x, y):
    return [b * y + x * a + e , d * y + x * c + f]


# In[57]:

"""
Align images and blend

a) Affine
""";


# In[58]:

a = float(z[0])
b = float(z[1])
c = float(z[2])
d = float(z[3])
e = float(z[4])
f = float(z[5])


# In[59]:

nlX = np.zeros_like(lX)
nlY = np.zeros_like(lY)


# In[60]:

# Calc new (affine transformed) points
for i in range(N): 
    nlX[i], nlY[i] = affine_transform_point(lX[i], lY[i])


# In[61]:

scatter = Image.new('RGB', (w, h), 'lightgray')


draw = ImageDraw.Draw(scatter)
draw.ellipse((x_c-ELL_RAD, y_c-ELL_RAD, x_c+ELL_RAD, y_c+ELL_RAD), fill='darkgreen')
# Draw central part boundary
draw.ellipse((x_c-CENTER_RAD, y_c-CENTER_RAD, x_c+CENTER_RAD, y_c+CENTER_RAD), outline='black')

for i in range(N): # draw star points
    draw.ellipse((nlX[i]-ELL_RAD, nlY[i]-ELL_RAD, nlX[i]+ELL_RAD, nlY[i]+ELL_RAD), fill='red')
    draw.ellipse((rX[i]-ELL_RAD, rY[i]-ELL_RAD, rX[i]+ELL_RAD, rY[i]+ELL_RAD), fill='blue')
# scatter


# In[ ]:




# In[ ]:




# In[62]:

"""
b) Affine + Ditortion3
""";


# In[63]:

a = float(z3[0])
b = float(z3[1])
c = float(z3[2])
d = float(z3[3])
e = float(z3[4])
f = float(z3[5])

eps1 = float(z3[6])
eps2 = float(z3[7])


# In[64]:

nlX3 = np.zeros_like(lX)
nlY3 = np.zeros_like(lY)
nrX3 = np.zeros_like(rX)
nrY3 = np.zeros_like(rY)


# In[65]:

# Get rid of distortions on new left Img
for i in range(N): 
    dist_l = (lX[i] - x_c)**2 + (lY[i] - y_c)**2

    zx1 = (lX[i] - x_c) * dist_l
    zy1 = (lY[i] - y_c) * dist_l


    nlX3[i] = lX[i] - eps1 * zx1
    nlY3[i] = lY[i] - eps1 * zy1


# Get rid of distortions on new right Img

for i in range(N): 
    dist_r = (rX[i] - x_c)**2 + (rY[i] - y_c)**2

    zx2 = (rX[i] - x_c) * dist_r
    zy2 = (rY[i] - y_c) * dist_r


    nrX3[i] = rX[i] - eps2 * zx2
    nrY3[i] = rY[i] - eps2 * zy2


# In[66]:

# Calc new (affine transformed) points
for i in range(N): 
    nlX3[i], nlY3[i] = affine_transform_point(lX[i], lY[i])


# In[67]:

scatter3 = Image.new('RGB', (w, h), 'lightgray')


draw = ImageDraw.Draw(scatter3)
draw.ellipse((x_c-ELL_RAD, y_c-ELL_RAD, x_c+ELL_RAD, y_c+ELL_RAD), fill='darkgreen')
# Draw central part boundary
draw.ellipse((x_c-CENTER_RAD, y_c-CENTER_RAD, x_c+CENTER_RAD, y_c+CENTER_RAD), outline='black')

for i in range(N): # draw star points
    draw.ellipse((nlX3[i]-ELL_RAD, nlY3[i]-ELL_RAD, nlX3[i]+ELL_RAD, nlY3[i]+ELL_RAD), fill='red')
    draw.ellipse((nrX3[i]-ELL_RAD, nrY3[i]-ELL_RAD, nrX3[i]+ELL_RAD, nrY3[i]+ELL_RAD), fill='blue')
# scatter3


# In[ ]:




# In[68]:

"""
c) Affine + Ditortion3,5
""";


# In[69]:

a = float(z35[0])
b = float(z35[1])
c = float(z35[2])
d = float(z35[3])
e = float(z35[4])
f = float(z35[5])

eps1 = float(z35[6])
eps2 = float(z35[7])
eps3 = float(z35[8])
eps4 = float(z35[9])


# In[70]:

nlX35 = np.zeros_like(lX)
nlY35 = np.zeros_like(lY)
nrX35 = np.zeros_like(rX)
nrY35 = np.zeros_like(rY)


# In[71]:

# Get rid of distortions on new left Img
for i in range(N): 
    dist_l = (lX[i] - x_c)**2 + (lY[i] - y_c)**2

    zx1 = (lX[i] - x_c) * dist_l
    zy1 = (lY[i] - y_c) * dist_l
    wx1 = (lX[i] - x_c) * dist_l**2
    wy1 = (lY[i] - y_c) * dist_l**2

    
    nlX35[i] = lX[i] - eps1 * zx1 - eps3 * wx1
    nlY35[i] = lY[i] - eps1 * zy1 - eps3 * wy1


# Get rid of distortions on new right Img

for i in range(N): 
    dist_r = (rX[i] - x_c)**2 + (rY[i] - y_c)**2

    zx2 = (rX[i] - x_c) * dist_r
    zy2 = (rY[i] - y_c) * dist_r
    wx2 = (rX[i] - x_c) * dist_r**2
    wy2 = (rY[i] - y_c) * dist_r**2


    nrX35[i] = rX[i] - eps2 * zx2 - eps4 * wx2
    nrY35[i] = rY[i] - eps2 * zy2 - eps4 * wy2


# In[72]:

# Calc new (affine transformed) points
for i in range(N): 
    nlX35[i], nlY35[i] = affine_transform_point(lX[i], lY[i])


# In[73]:

scatter35 = Image.new('RGB', (w, h), 'lightgray')


draw = ImageDraw.Draw(scatter35)
draw.ellipse((x_c-ELL_RAD, y_c-ELL_RAD, x_c+ELL_RAD, y_c+ELL_RAD), fill='darkgreen')
# Draw central part boundary
draw.ellipse((x_c-CENTER_RAD, y_c-CENTER_RAD, x_c+CENTER_RAD, y_c+CENTER_RAD), outline='black')

for i in range(N): # draw star points
    draw.ellipse((nlX35[i]-ELL_RAD, nlY35[i]-ELL_RAD, nlX35[i]+ELL_RAD, nlY35[i]+ELL_RAD), fill='red')
    draw.ellipse((nrX35[i]-ELL_RAD, nrY35[i]-ELL_RAD, nrX35[i]+ELL_RAD, nrY35[i]+ELL_RAD), fill='blue')
# scatter35


# In[ ]:




# Выходные данные программы:
# 
# 1) $\Delta x_i, \Delta y_i, \; i = 1,N$
# 
# 2) $\sigma^2 = \frac{1}{N} \sum\limits_{i=1}^{N} 
#                 \left( \Delta x_i^2 + \Delta y_i^2 \right)$
#                 
# 3) $ 
#      M_x = \max\limits_{i=1,N} \{\Delta x_i\}, \;
#      M_y = \max\limits_{i=1,N} \{\Delta y_i\}, \;
#      M = \max\{M_x, M_y\}
#    $                

# In[74]:

"""
For all stars -- calculate distance (by x, by y) 
between corresponding stars on aligned images:
\Delta_x, \Delta_y;
\sigma^2;
M_x, M_y, M
""";


# In[75]:

"""
a) Affine
""";


# In[76]:

delX = np.zeros_like(nlX)
delY = np.zeros_like(nlY)


# In[78]:

delX = abs(nlX - rX)
delY = abs(nlY - rY)
print("delX:", delX)
print("delY:", delY)


# In[79]:

sigSqr = 1.0 / N * sum(delX**2 + delY**2)
mX = max(delX)
mY = max(delY)
m = max(mX, mY)

print("mX: %.4f mY: %.4f m: %.4f" % (mX, mY, m))
print("sigSqr: %.4f" % sigSqr)


# In[80]:

"""
b) Affine + Ditortion3
""";


# In[81]:

delX3 = np.zeros_like(nlX3)
delY3 = np.zeros_like(nlY3)


# In[82]:

delX3 = abs(nlX3 - nrX3)
delY3 = abs(nlY3 - nrY3)
print("delX3:", delX3)
print("delY3:", delY3)


# In[83]:

sigSqr3 = 1.0 / N * sum(delX3**2 + delY3**2)
mX3 = max(delX3)
mY3 = max(delY3)
m3 = max(mX3, mY3)
print("mX3: %.4f mY3: %.4f m3: %.4f" % (mX3, mY3, m3))

print("sigSqr3: %.4f" % sigSqr3)


# In[84]:

"""
b) Affine + Ditortion3,5
""";


# In[85]:

delX35 = np.zeros_like(nlX35)
delY35 = np.zeros_like(nlY35)


# In[86]:

delX35 = abs(nlX35 - nrX35)
delY35 = abs(nlY35 - nrY35)
print("delX35:", delX35)
print("delY35:", delY35)


# In[87]:

sigSqr35 = 1.0 / N * sum(delX35**2 + delY35**2)
mX35 = max(delX35)
mY35 = max(delY35)
m35 = max(mX35, mY35)

print("mX35: %.4f mY35: %.4f m35: %.4f" % (mX35, mY35, m35))
print("sigSqr35: %.4f" % sigSqr35)


# In[ ]:




# In[88]:

"""
Save results: mode=("a" / "d3" / "d35")
radius, N, sigma^2, M
""";


# In[89]:

res_fname = results_folder + os.sep + "results_a.txt"
res_fname3 = results_folder + os.sep + "results_d3.txt"
res_fname35 = results_folder + os.sep + "results_d35.txt"


# In[90]:

print(res_fname)
print(res_fname3)
print(res_fname35)


# In[91]:

header = "# CENTER_RAD, N, sigSqr, m"
line = " ".join(str(e) for e in [CENTER_RAD, N, sigSqr, m])
print(line)


# In[92]:

line3 = " ".join(str(e) for e in [CENTER_RAD, N, sigSqr3, m3])
print(line3)


# In[93]:

line35 = " ".join(str(e) for e in [CENTER_RAD, N, sigSqr35, m35])
print(line35)


# In[94]:

def write_to_file(res_fname, line):
    # if file exists -- append
    if os.path.isfile(res_fname):
        with open(res_fname, 'a') as fout:
            fout.write(line)
            fout.write('\n')

    # If file doest exist -- create new
    else:
        with open(res_fname, 'w') as fout:
            fout.write(header)
            fout.write('\n')
            fout.write(line)
            fout.write('\n')


# In[95]:

write_to_file(res_fname, line)
write_to_file(res_fname3, line3)
write_to_file(res_fname35, line35)


# In[ ]:




# "4." Построить графики зависимости $\sigma$ и $M$ от размера круга для аффинного, аффин+дист.3 порядка и аффин+дист.3 порядка+дист 5
# порядка.

# In[96]:

"""
Plot graph sigma^2(radius), M(radius)
""";


# In[97]:

# Load data from files "results_a.txt, results_d3.txt, results_d35.txt"
# (without header -- 1st string)
# load to arrays resA, resD3, resD35


# In[98]:

# Make data unique in resA, resD3, resD35 
# (remove duplicating rows)


# In[99]:

# Make arrays of radius, Num stars, sigma Squared, M


# In[100]:

# Plot sigma^2(radius), M(radius) with 3 curves on each (A, D3, D35)


# In[ ]:



