# This utility copies JS plugins from the repository to the Reaper JS plugin folder
# This is handy when developing from within the local repository in order to test in Reaper


require 'FileUtils'


# Determine what platform we are on
def mac?
  (Object::RUBY_PLATFORM =~ /darwin/i) ? true : false
end

def linux?
   (Object::RUBY_PLATFORM =~ /linux/i) ? true : false
end

def win?
  !mac? && !linux?
end



# Copy plugins to the appropriate Reaper folder
def copyPluginsToReaper
  
  if mac?
    
    # Create plugin folder if it doesn't exist already
    pluginFolderPath = File.expand_path("~") + "/Library/Application Support/REAPER/Effects/AmbiToolkit"
    FileUtils.mkdir_p(pluginFolderPath) unless File.exists?(pluginFolderPath)
    
    # Copy plugins
    `cp -r "../plugins/" "#{pluginFolderPath}"`
  end
  
  if win?
    puts "TODO: Need to add support for Windows in this script"
  end
  
end


# Run the method
copyPluginsToReaper
