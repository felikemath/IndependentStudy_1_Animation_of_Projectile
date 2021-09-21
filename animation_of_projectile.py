from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib.patches import Arc
import math


x0 = 0.0
y0 = 0.0
v0 = 20  # m/s
theta0 = 30  # initial launching angle in degrees
k = 0.1  # linear drag force coefficient of air resistance

# gravity constant g
g = 9.8 # m/s^2
m = 1.0 # kg
# derived variables for convenience
cos_theta = math.cos(theta0 * math.pi / 180.0)
sin_theta = math.sin(theta0 * math.pi / 180.0)
vx0 = v0 * cos_theta
vy0 = v0 * sin_theta

# A: without air resistance; B: with air resistance
prev_xA = 0
prev_yA = 0
prev_distA = 0
prev_xB = 0
prev_yB = 0
prev_distB = 0

xArrayA = []
xArrayB = []
yArrayA = []
yArrayB = []

intervalms = 20  # this means 20 ms per frame
numframes = 300
duration = intervalms*numframes/1000.0  # sec

# calculate the highest vertical point & lowest vertical points during the duration of (numframes * intervalms)
highestpoint = vy0*vy0/(2.0*g)+y0
lowestpoint = vy0 * duration - 0.5*g*duration*duration+y0
if lowestpoint > y0:
    lowestpoint = y0
scale_arrow = (highestpoint-lowestpoint) / v0 * 0.1
x_range = vx0*duration
y_range = highestpoint+5 - lowestpoint
if y_range > x_range:
    x_range = y_range

# Set up the figure, axis, and the plot elements to animate
fig = plt.figure(figsize=(8, 8))
ax = plt.axes(xlim=(x0, x_range+x0), ylim=(lowestpoint, highestpoint+5), xlabel='x position (m)', ylabel='y position (m)')

plt.title("Projectile Motion with and without Air Resistance")
ax.arrow(x0, y0, scale_arrow*vx0, scale_arrow*vy0, head_width=1, head_length=2, fc='g', ec='g')
if theta0 >= 30:
    plt.text(x0+x_range/20, y0+x_range/20*math.sin(theta0*math.pi/360), r'$\theta_0 $', color='black', fontsize=14, ha='left', va='center')
else:
    plt.annotate(r'${\theta_0}$', xy=(5, y0+x_range/20*math.sin(theta0*math.pi/360)), xycoords='data',
                 xytext=(x0+x_range/20, y0-y_range/20), textcoords='data', fontsize=14, color='black',
                 arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))

ax.plot([x0, x_range/10], [y0, y0], '--k', lw=0.5)

dotB, = ax.plot([], [], 'or', lw=0.4)
lineB, = ax.plot([], [], '--m', lw=1, label='w/ air resistance')
dotA, = ax.plot([], [], 'ob', lw=0.4)
lineA, = ax.plot([], [], '--b', lw=1, label='w/o air resistance')

# The following lines help find appropriate empty locations to display the annotations
if theta0 >= 60:
    x_display_loc = 0.5
    y_display_loc = 0.8
else:
    x_display_loc = 0.02
    y_display_loc = 0.5
y_display_loc_inc = 0.03

text_drag_k = ax.text(x_display_loc, y_display_loc + 1*y_display_loc_inc, '', transform=ax.transAxes)
text_init_x = ax.text(x_display_loc, y_display_loc - 0*y_display_loc_inc, '', transform=ax.transAxes)
text_init_y = ax.text(x_display_loc, y_display_loc - 1*y_display_loc_inc, '', transform=ax.transAxes)
text_init_angle = ax.text(x_display_loc, y_display_loc - 2*y_display_loc_inc, '', transform=ax.transAxes)
text_init_speed = ax.text(x_display_loc, y_display_loc - 3*y_display_loc_inc, '', transform=ax.transAxes)
text_time = ax.text(x_display_loc, y_display_loc - 5*y_display_loc_inc, '', transform=ax.transAxes)
text_x = ax.text(x_display_loc, y_display_loc - 6*y_display_loc_inc, '', transform=ax.transAxes)
text_y = ax.text(x_display_loc, y_display_loc - 7*y_display_loc_inc, '', transform=ax.transAxes)
text_vx = ax.text(x_display_loc, y_display_loc - 8*y_display_loc_inc, '', transform=ax.transAxes)
text_vy = ax.text(x_display_loc, y_display_loc - 9*y_display_loc_inc, '', transform=ax.transAxes)
text_speed = ax.text(x_display_loc, y_display_loc - 10*y_display_loc_inc, '', transform=ax.transAxes)
text_dist = ax.text(x_display_loc, y_display_loc - 11*y_display_loc_inc, '', transform=ax.transAxes)
text_energy_k = ax.text(x_display_loc, y_display_loc - 12*y_display_loc_inc, '', transform=ax.transAxes)
text_energy_p = ax.text(x_display_loc, y_display_loc - 13*y_display_loc_inc, '', transform=ax.transAxes)
text_energy_m = ax.text(x_display_loc, y_display_loc - 14*y_display_loc_inc, '', transform=ax.transAxes)

ax.set_aspect('equal', adjustable='box')
leg = ax.legend()
ax.add_patch(Arc((0, 0), 5, 5, theta1=0.0, theta2=theta0, edgecolor='black', lw=1.5))
# fig.tight_layout()


# initialization function: automatically called once before the 1st frame
def init():
    global prev_distA, prev_distB, prev_xA, prev_yA, prev_xB, prev_yB
    prev_xA = x0
    prev_yA = y0
    prev_xB = x0
    prev_yB = y0
    prev_distA = 0
    prev_distB = 0

    xArrayA[:] = []
    yArrayA[:] = []
    xArrayB[:] = []
    yArrayB[:] = []
    lineA.set_data([], [])
    dotA.set_data([], [])
    lineB.set_data([], [])
    dotB.set_data([], [])

    text_time.set_text('')
    text_x.set_text('')
    text_y.set_text('')
    text_vx.set_text('')
    text_vy.set_text('')
    text_speed.set_text('')
    text_dist.set_text('')
    text_energy_k.set_text('')
    text_energy_p.set_text('')
    text_energy_m.set_text('')

    text_drag_k.set_text('Linear air resistance k = {:.3f} /s'.format(k))
    text_init_x.set_text('initial x loc = {:.1f} m'.format(x0))
    text_init_y.set_text('initial y loc = {:.1f} m'.format(y0))
    text_init_angle.set_text(r'initial angle $\theta_0  = {:.1f}^\circ $'.format(theta0))
    text_init_speed.set_text('initial speed = {:.1f} m/s'.format(v0))

    return dotA, lineA, dotB, lineB, text_time, text_y, text_x, text_vy, text_vx, text_speed, text_dist, text_energy_k, text_energy_p, text_energy_m


# animation function
def animate(i):
    global prev_distA, prev_distB, prev_xA, prev_yA, prev_xB, prev_yB

    # time sampling
    t = i * intervalms / 1000.0

    # A: without air resistance
    xA = x0 + vx0 * t
    yA = y0 + vy0 * t - (g*t*t)/2.0
    vxA = vx0
    vyA = vy0 - g*t
    vA = math.sqrt(vxA * vxA + vyA * vyA)

    # B: with air resistance
    xB = x0 + vx0 * (1 - math.exp(-k*t/m)) * m / k
    yB = y0 + vy0 * (1 - math.exp(-k*t/m)) * m / k + g * (1 - k*t/m - math.exp(-k*t/m)) * (m*m)/(k*k)
    vxB = vx0 * math.exp(-k*t/m)
    vyB = vy0 * math.exp(-k*t/m) - g * (1 - math.exp(-k*t/m)) * m/k
    vB = math.sqrt(vxB * vxB + vyB * vyB)

    xArrayA.append(xA)
    yArrayA.append(yA)
    lineA.set_data(xArrayA, yArrayA)
    dotA.set_data(xA, yA)

    xArrayB.append(xB)
    yArrayB.append(yB)
    lineB.set_data(xArrayB, yArrayB)
    dotB.set_data(xB, yB)

    text_time.set_text('time = %.2f s' % t)
    text_x.set_text('x ($m$):      w/ = {:.3f},  w/o = {:.3f}'.format(xB, xA))
    text_y.set_text('y ($m$):      w/ = {:.3f},  w/o = {:.3f}'.format(yB, yA))

    text_vx.set_text('Vx ($m/s$): w/ = {:.3f},  w/o = {:.3f}'.format(vxB, vxA))
    text_vy.set_text('Vy ($m/s$): w/ = {:.3f},  w/o = {:.3f}'.format(vyB, vyA))

    text_speed.set_text('Speed ($m/s$): w/ = {:.3f}, w/o = {:.3f}'.format(vB, vA))

    distA = prev_distA + math.sqrt((xA-prev_xA)**2+(yA-prev_yA)**2)
    distB = prev_distB + math.sqrt((xB - prev_xB) ** 2 + (yB - prev_yB) ** 2)

    text_dist.set_text('TravelDistance ($m$): w/ = {:.1f}, w/o = {:.1f}'.format(distB, distA))

    prev_xA = xA
    prev_yA = yA
    prev_xB = xB
    prev_yB = yB
    prev_distA = distA
    prev_distB = distB

    e_kinetic_A = 0.5 * m * vA * vA
    e_kinetic_B = 0.5 * m * vB * vB
    e_potential_A = m * g * yA
    e_potential_B = m * g * yB
    eA = e_kinetic_A + e_potential_A
    eB = e_kinetic_B + e_potential_B
    text_energy_k.set_text('K.Energy ($J$): w/ = {:.1f}, w/o = {:.1f}'.format(e_kinetic_B, e_kinetic_A))
    text_energy_p.set_text('P.Energy ($J$): w/ = {:.1f}, w/o = {:.1f}'.format(e_potential_B, e_potential_A))
    text_energy_m.set_text('M.Energy ($J$): w/ = {:.1f}, w/o = {:.1f}'.format(eB, eA))

    return dotA, lineA, dotB, lineB, text_time, text_y, text_x, text_vy, text_vx, text_speed, text_dist, text_energy_k, text_energy_p, text_energy_m


def main():
    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=numframes, interval=intervalms, blit=True)
    # save the animation as an mp4. Please install ffmpeg and add its bin folder to PATH
    anim.save(r'output\Projectile_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
    plt.show()


if __name__=='__main__':
    main()