# plotting ambisonic transforms
# set sizes with actual dimensions
# include layouts in appropriate scalings

# debug. . .

# matplotlib, muse

import numpy as np
import muse as mu
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import colors, colorbar


# included here as there is an error of some sort with analysis.py
# NOT SURE if this is as a result of moving to EPD
# vaed (azimuth, elevation, directivity)
def vaed(b):
    """vaed(a)
    
    Analyze an ambisonic B-format sound field, returning azimuth,
    elevation and directivity. Returns values as time varying.
    
    Inputs:
        - b         : Input b-format signal

    Outputs: ([a, e, d])
    
      [a, e, d] -- Azimuth, elevation and directivity in radians.
                   (See direct for details on directivity.)

    """

    # calculate directivity
    b_sqrd = b**2

    d = 2 * mu.arctan2(
        mu.rec_sqrt2 * mu.sqrt(
            (b_sqrd[:, 1] + b_sqrd[:, 2] + b_sqrd[:, 3])
            ),
        mu.sqrt(b_sqrd[:, 0])
        )

##    # adjust directivity so that it is idealized
##    b_dir = mu.direct(
##        b,
##        2. * mu.arctan(1 / mu.tan(d / 2.))
##        )
##
##    # calculate azimuth, elevation
##    b_abs = mu.a_to_b(
##        mu.abs(
##            mu.b_to_a(b_dir, weight = 'car')
##            ),
##        weight = 'car'
##        )

    # calculate azimuth, elevation
    b_abs = mu.a_to_b(
        mu.abs(
            mu.b_to_a(b, weight = 'car')
            ),
        weight = 'car'
        )


    # translated to [r, theta, phi]
    spher = mu.cart_to_spher(b_abs[:, 1:])

    # return [theta, phi, d]
    res = mu.hstack(
        (spher[:, 1:],
        mu.interleave(d))
        )

    return res


# calculate gain
# not quite sure if this gives energy or velocity
def calc_gain(b):

    # calculate amp
    res = np.sqrt(
        np.sum(
            np.square(
                np.array([mu.sqrt2, 1., 1., 1.]) * b
                ), axis = 1
            ) / 2.
        )

    # return gain
    return mu.amp_to_db(res)


# linear rescale
def rescale(x, x_min_max = (0., 1.), y_min_max = (0., 1.)):

    m = (y_min_max[1] - y_min_max[0]) / (x_min_max[1] - x_min_max[0])
    res = y_min_max[0] + (x - x_min_max[0]) * m

    return res

# cm / inch conversions
def in_to_cm(x):
    return 2.54 * x

def cm_to_in(x):
    return x/in_to_cm(1.)


# ********************************************************
# set up params:

# display?
fig_display = False
# fig_display = True

# save image?
fig_save = True
# fig_save = False
# fig_format = 'pdf'
# fig_format = 'svg'
fig_format = 'png'
fig_dir = '/Users/josephla/Desktop/new ATK pics/'
# fig_dir = '/Users/josephla/Desktop/'
# fig_dir = '/Users/josephla/Desktop/SBCM paper/tex/sbc-JA-latex2003/'
# fig_tparent = True              # transparent?
fig_tparent = False              # transparent?

trans_list = [
    (mu.focus_x, 'Focus', 'focus_fig'),
    (mu.push_x, 'Push', 'push_fig'),
    (mu.press_x, 'Press', 'press_fig'),
    (mu.zoom_x, 'Zoom', 'zoom_fig'),
    (mu.zoom_y, 'Balance', 'balance_fig'),
    (mu.asymmetry_x, 'Asymmetry', 'asymmetry_fig'),
    (mu.direct, 'Direct', 'direct_fig'),
    (mu.direct_x, 'DirectY', 'direct_x_fig'),
    (mu.direct_y, 'DirectY', 'direct_y_fig')
]

# select transform from above list:
# t_indx = 0
# t_indx = 1
# t_indx = 2
# t_indx = 3
# t_indx = 4
t_indx = 5
# t_indx = 6
# t_indx = 7
# t_indx = 8

transfunc, suptitle, fig_file = trans_list[t_indx]

# plot colors, prettines. . .
cmap = cm.jet
# cmap = cm.binary
# grid_color = '#2937f2'
grid_color = 'gray'
lin_color = '#ee8d18'

# dpi = 300
# dpi = 100       # 1190 x 313 pixels
# dpi = 60        # 714 x 188 pixels
# dpi = 54        # 642 x 169 pixels
dpi = 56        # 666 x 175 pixels


# params to generate plot data. . .
lin_points = 360                # points to 'scatter plot' for 'continuous' image
gain_range = (-24., 6.)         # min/max vals to map gain against
db_min = -180.                  # min value to plot against (used for dominance)
c_points_num = 8                # number of cardinal points to plot

t_angle_num = 4                 # number of transform angles to plot (=num plots)
# t_angle_num = 5                 # number of transform angles to plot (=num plots)
# t_angle_num = 7                 # number of transform angles to plot (=num plots)

t_angle_min = 0.                # min angle in radians
# t_angle_max = mu.halfPi         # max "     "  "
t_angle_max = -mu.halfPi         # max "     "  "
t_angles = mu.lin([t_angle_min, t_angle_max], t_angle_num)


# *************
# parms for setting up page. . .

# dimensions in cm
plot_size = cm_to_in(6.)       # width & height
wtop = cm_to_in(1.)             # white space: top
wbottom = wtop                  # white space: bottom
# wsides = cm_to_in(1.)           # white space: left/right side
wsides = cm_to_in(0.75)           # white space: left/right side
wspace = cm_to_in(.75)           # white space: between subplots
cbar_width = cm_to_in(2.5)        # width

# set up necessary scalings. . .
fig_height = wtop + wbottom + plot_size
fig_width = (2 * wsides) + (t_angle_num * plot_size) + ((t_angle_num - 1) * wspace) + cbar_width


print fig_width, fig_height
print in_to_cm(fig_width), in_to_cm(fig_height)

# *************
# parms for setting up plotting formats. . .
# set grid lines
plt.rc('grid', color = grid_color, linewidth=1, linestyle='--')

# figsize is in (width, height) in inches
fig = plt.figure(figsize = (fig_width, fig_height))

# set up subplot size, position. . .
# can use left and right to scale the images between left and right boundaries of image
# fig.subplots_adjust(left=None, bottom=None, right=None, wspace=None, hspace=None)
# fig.SubplotParams(left=.125, right=.9, bottom=.1, top = .9, wspace=.2, hspace=.2)
fig.subplots_adjust(
    left = (wsides + cbar_width) / fig_width,
    right = 1. - (wsides / fig_width),
    bottom = wbottom / fig_height,
    top = 1. - (wtop / fig_height),
    wspace = wspace
)


# titles & labels
axis_labels = ('', '', '+X','', '+Y', '', '', '')

# TITLE quoted out. . .
# suptitle(t, **kwargs) Add a centered title to the figure
# fig.suptitle(suptitle, fontsize=30, fontweight = 'regular')

add_title = False
if add_title:
    fig.suptitle(suptitle, fontsize=30, fontweight = 'regular')


# misc ploting data
plt_angle = mu.halfPi           # use to rotate image for plotting: X on vertical


# annotations of selected cardinal points:
# Front Centre, Back Centre, Left, Back Left

# some defaults:
textcoords = 'axes fraction'    # use axes fraction when these are separate axes
arrowstyle = "simple"
facecolor = "0.3"
edgecolor = "none"
# shrinkB = 5.
shrinkB = 10.
an_font = 'sans-serif'
# an_fontsize = 8
an_fontsize = 10

# params
display_gain = True
display_elev = False


# ***********************************************************
# set up and generate colorbar
# includes initial set up for subplotting

# cmap = cm.jet               # good: close to audio metering
norm = colors.Normalize(gain_range[0], gain_range[1], clip = True)

cax = fig.add_axes([
        wsides/fig_width + (.33 * cbar_width/fig_width),
        wbottom/fig_height + (.1 * plot_size / fig_height),
        .125 * cbar_width/fig_width,
        .8 * plot_size/fig_height
])


cb1 = colorbar.ColorbarBase(
    cax,
    cmap = cmap,
    norm = norm
    )

# cb1.set_label('gain (dB)', fontsize = 8)
cb1.set_label('gain (dB)', fontsize = 10)

# cb1.ax is matplotlib.axes.Axes object
# get the yaxis of the colorbar
# move ticks and label to left
# do expect there is a way to do so when settin gup ColorbarBase

cb1.ax.get_yaxis().set_ticks_position('left')
cb1.ax.get_yaxis().set_label_position('left')
# below is a hack. . . better to get it working automatically, but how?
cb1.ax.get_yaxis().set_ticklabels(
    ['-24', '', '-18', '', '-12', '', '-6', '', '0', '', '+6'],
#     fontsize = 8
    fontsize = 10
    )


# ***********************************************************
# generate data to plot

# the idea is this:
# 1) generate a 'unity' signal in b-format: a_bf
# 2) transform that signal via ATK transform: b_bf
# 3) return r & theta from transformed signal: 
# 4) plot!!

# b-format test 'signal'
# the 'signal' is actually merely an envelope (of 1s)
a_theta =  mu.twoPi * np.arange(lin_points) / lin_points
a_bf = mu.mono_to_b(mu.ones(lin_points), a_theta)

# empty axes list: populate with subplots
axs = []

# ***********************************************************
# big loop for multi-plots

for (subplot_num, t_angle) in enumerate(t_angles):

    print "plotting subplot ", subplot_num

    # ***********************************************************
    # transform the b-format signal
    b_bf = apply(transfunc, (a_bf, t_angle))

    # calculate values for plotting

    #  gain of transformed b-format
    b_gain = calc_gain(b_bf)

    # calculate theta using vaed
    # should use mu.vaed, however, there is an error to be corrected in muse
    b_vaed = vaed(b_bf) # calulate azimuth, elevation, directivity
    b_theta = b_vaed[:,0]

    # calculate r projected in the XY plane
    # normalized by W: looks good, use this!
    b_r = mu.deinterleave(mu.cart_to_pol(b_bf[:,1:3] / (mu.sqrt2 * mu.interleave(b_bf[:,0]))))[0]

    # calculate Z normalized by W
    b_z = b_bf[:,3] / (mu.sqrt2 * b_bf[:,0])

    # calculate elevation, phi
    # normalized by W
    b_elev = mu.deinterleave(mu.cart_to_spher(b_bf[:,1:] / (mu.sqrt2 * mu.interleave(b_bf[:,0]))))[2]

    # generate data for cardinal points from b_theta, b_r, b_gain
    c_indx = np.arange(0, lin_points, lin_points/c_points_num)

    c_gain = np.empty((c_points_num), float)
    c_theta = np.empty((c_points_num), float)
    c_r = np.empty((c_points_num), float)
    c_z = np.empty((c_points_num), float)
    c_elev = np.empty((c_points_num), float)

    # fill c_theta, c_r, c_gain
    for (i, indx) in enumerate(c_indx):
        c_gain[i] = b_gain[indx]
        c_theta[i] = b_theta[indx]
        c_r[i] = b_r[indx]
        c_z[i] = b_z[indx]
        c_elev[i] = b_elev[indx]


    # ***********************************************************
    # generate the plots

    # add subplot
    # subplot(numRows, numCols, plotNum)
    axs.append(fig.add_subplot(1, t_angle_num, (subplot_num + 1), projection = 'polar'))

    # plot angular distortion: as a scatter plot
    axs[subplot_num].scatter(
        plt_angle + b_theta,
        b_r,
        s = 1,
        color = cmap(rescale(b_gain, gain_range)),
        alpha = 1.
        )

    # plot cardinal points: as scatter plot
    axs[subplot_num].scatter(
        plt_angle + c_theta,
        c_r,
#         s = rescale(c_gain, gain_range, (1., 100.)),
        s = rescale(c_gain, gain_range, (1., 75.)),
        color = cmap(rescale(c_gain, gain_range)),
#         alpha = .75
        alpha = 1.
        )


# ************************************************************
# annotations of selected cardinal points:
# Front Centre, Back Centre, Left, Back Left

    # set labels
    if display_gain is True:
        fc_text = 'Front\n%.1f dB' % c_gain[0]
        bc_text = 'Back\n%.1f dB' % c_gain[4]
        sl_text = 'Left\n%.1f dB' % c_gain[2]
        bl_text = 'Back-Left\n%.1f dB' % c_gain[3]
#         fc_text = '%.1f dB' % c_gain[0]
#         bc_text = '%.1f dB' % c_gain[4]
#         sl_text = '%.1f dB' % c_gain[2]
#         bl_text = '%.1f dB' % c_gain[3]

    elif display_elev is True:
        fc_text = 'Front\nz = %.1f' % c_z[0]
        bc_text = 'Back\nz = %.1f' % c_z[4]
        sl_text = 'Left\nz = %.1f' % c_z[2]
        bl_text = 'Back-Left\nz = %.1f' % c_z[3]
#         fc_text = 'z = %.1f' % c_z[0]
#         bc_text = 'z = %.1f' % c_z[4]
#         sl_text = 'z = %.1f' % c_z[2]
#         bl_text = 'z = %.1f' % c_z[3]

    else:
        fc_text = 'Front'
        bc_text = 'Back'
        sl_text = 'Left'
        bl_text = 'Back-Left'

    # annotate position of Front Centre
    axs[subplot_num].annotate(
        fc_text,
        xy = (plt_angle + c_theta[0], c_r[0]),  # theta, radius
#         xytext = (1., 1.),    # text position as fraction, fraction
        xytext = (1.1, 1.05),    # text position as fraction, fraction
        textcoords = textcoords,
        arrowprops = dict(
            arrowstyle = arrowstyle,
            facecolor = facecolor,
            edgecolor = edgecolor,
            connectionstyle = ((
                    "arc3,rad=0.35",
                    "arc3,rad=-0.35")[c_theta[0] < -mu.halfPi/2]),
                       shrinkB = shrinkB),
        horizontalalignment='right',
        verticalalignment='bottom',
        fontname = an_font,
        fontsize = an_fontsize
        )

    # annotate position of Back Centre
    if c_gain[4] > db_min:       # check if gain is greater than -00
        axs[subplot_num].annotate(
            bc_text,
            xy = (plt_angle + c_theta[4], c_r[4]),  # theta, radius
#             xytext = (1., 0.),    # text position as fraction, fraction
            xytext = (1.1, -.1),    # text position as fraction, fraction
            textcoords = textcoords,
            arrowprops = dict(
                arrowstyle = arrowstyle,
                facecolor = facecolor,
                edgecolor = edgecolor,
                connectionstyle = (
#                     (
#                         "arc3,rad=-0.35",
#                         "arc3,rad=0.35")[(
#                             (c_theta[4] >= mu.pi + mu.halfPi/2) or
#                             (c_theta[4] <= .01) or
#                             (c_r[4] <= .01))]),
#                         "arc3,rad=0.35"),
                        "arc3,rad=-0.35"),
                shrinkB = shrinkB),
            horizontalalignment='right',
            verticalalignment='bottom',
            fontname = an_font,
            fontsize = an_fontsize
            )
    else:
        axs[subplot_num].annotate(
#             'Back Centre\n-inf dB',
#             '-inf dB',
            'Back\n-inf dB',
            xy = (plt_angle + c_theta[4], c_r[4]),  # theta, radius
#             xytext = (1., 0.),    # text position as fraction, fraction
            xytext = (1.1, -.1),    # text position as fraction, fraction
            textcoords = textcoords,
            arrowprops = dict( # invisible
                facecolor = "none",
                edgecolor = "none",
                ),
            horizontalalignment='right',
            verticalalignment='bottom',
            fontname = an_font,
            fontsize = an_fontsize
            )        

    # annotate position of Side Left
    axs[subplot_num].annotate(
        sl_text,
        xy = (plt_angle + c_theta[2], c_r[2]),  # theta, radius
#         xytext = (0., 1.),    # text position as fraction, fraction
        xytext = (-.1, 1.05),    # text position as fraction, fraction
        textcoords = textcoords,
        arrowprops = dict(
            arrowstyle = arrowstyle,
            facecolor = facecolor,
            edgecolor = edgecolor,
            connectionstyle = ((
                    "arc3,rad=0.35",
                    "arc3,rad=-0.35")[c_theta[2] < mu.halfPi/2]),
            shrinkB = shrinkB),
        horizontalalignment='left',
        verticalalignment='bottom',
        fontname = an_font,
        fontsize = an_fontsize
        )

    # annotate position of Back Left
    axs[subplot_num].annotate(
        bl_text,
        xy = (plt_angle + c_theta[3], c_r[3]),  # theta, radius
#         xytext = (0., 0.),    # text position as fraction, fraction
        xytext = (-.1, -.1),    # text position as fraction, fraction
        textcoords = textcoords,
        arrowprops = dict(
            arrowstyle = arrowstyle,
            facecolor = facecolor,
            edgecolor = edgecolor,
            connectionstyle = ((
                    "arc3,rad=0.35",
                    "arc3,rad=-0.35")[c_theta[3] <= mu.pi - mu.halfPi/2]),
            shrinkB = shrinkB),
        horizontalalignment='left',
        verticalalignment='bottom',
        fontname = an_font,
        fontsize = an_fontsize
        )

    # ************************************************************
    # then. . .
    # set up plotting axis range: limits of display
    # ylim appears to be the relevant value to set for polar
    axs[subplot_num].set_ylim(0., 1.2)

    # set up thetagrids
    # set the locations and labels of the radial gridlines and labels
    # might be able to move this outside the loop, just below
    plt.thetagrids(range(0,360,45), axis_labels,\
#                        frac = 1.2, fontname = 'serif', fontstyle = 'italic', fontsize = 12)
                       frac = 1.2, fontname = 'serif', fontstyle = 'italic', fontsize = 16)
#     plt.rgrids(mu.lin([0., 1.], 3)[1:], angle = 225., va = 'bottom', fontname = 'sanserif',\
#                    fontsize = 8)
    plt.rgrids(mu.lin([0., 1.], 3)[1:], angle = 247.5, va = 'bottom', fontname = 'sanserif',\
                   fontsize = 10)

    # then, set up title: may need to move this. . .
#     fig_title = '$\\omega = ' + str(mu.rad_to_deg(t_angle))[:-2] + '\degree$' # include degree symbol
    fig_title = '$\\alpha = ' + str(mu.rad_to_deg(t_angle))[:-2] + '\degree$' # include degree symbol
#     axs[subplot_num].set_title(fig_title, fontsize = 14,\
    axs[subplot_num].set_title(fig_title, fontsize = 18,\
                                   position = (.5, -.25)) # might want to specify this a bit better


# ************************************************
# display
if fig_display:
    plt.show()

# saving/exporting. . .
if fig_save:
	plt.savefig(fig_dir + fig_file + '.' + fig_format, format = fig_format, \
		transparent = fig_tparent, dpi = dpi)
