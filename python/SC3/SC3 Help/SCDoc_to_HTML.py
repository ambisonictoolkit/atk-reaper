# ****************************************************************************
# Massage help files generted by SCDoc for 3.4 release:
# 
# remove 'filename' element
# remove 'doclink' element
# remove 'inheritance' element (from)
# remove 'inheritance' element (to)
# remove 'inherited class methods'
# remove 'inherited instance methods'
# remove 'old help'
#
# ****************************************************************************

import os
import shutil
import lxml.html
import lxml.etree

# dirs
muse_dir = '/Users/josephla/Documents/Developer/ATK-Muse'

read_dir = '/sc/SC3ATK-help/SCDocHTMLHelp'
write_dir = '/sc/SC3ATK-help/Help'


# Help dirs
# NOTE: we could use os.walk, but this is slightly easier
help_dirs = os.listdir(muse_dir + read_dir)

# create write help dirs if they don't exist
for help_dir in help_dirs:
    if not os.path.isdir(os.path.join(muse_dir + write_dir, help_dir)):
        os.mkdir(os.path.join(muse_dir + write_dir, help_dir))


# itterate through files!
for help_dir in help_dirs:
    for file_name in os.listdir(os.path.join(muse_dir + read_dir, help_dir)):

        # read and write paths
        read_path = os.path.join(muse_dir + read_dir, help_dir, file_name)
        write_path = os.path.join(muse_dir + write_dir, help_dir, file_name)

        if os.path.splitext(read_path)[1] != '.html':

            # if not HTML, just copy
            shutil.copy(read_path, write_path)

        else:

            # else, process the HTML, using lxml

            # load the html
            htmltree = lxml.html.parse(read_path)

            # get the root
            root = htmltree.getroot()


            # NOTE: all the below code can likely be cleaned up
            #       with a better understanding of lxml!

            try: 
                # remove 'filename' element
                element = root.get_element_by_id('filename')
                element.getparent().remove(element)
            except:
                print 'No filename element!'

            # remove:
            #   'doclink'
            #   'inheritance' element (from)
            #   'inheritance' element (to)
            for element_class in ['doclink', 'inheritance', 'inheritance']:
                el_class = root.find_class(element_class)
                if len(el_class) != 0:
                    element = el_class[0]
                    element.getparent().remove(element)

            # remove:
            #   'inherited class methods'
            #   'inherited instance methods'
            #   'old help'
            for element in root.iter('*'):
                    if element.text == 'Inherited class methods' or \
                       element.text == 'Inherited instance methods' or \
                       element.text == 'old help':
                            element.getparent().remove(element)

            # write processed HTML out!
            write_file = open(write_path, 'w')
            write_file.write(
                lxml.etree.tostring(root, method='html', pretty_print=True))
            write_file.close()
