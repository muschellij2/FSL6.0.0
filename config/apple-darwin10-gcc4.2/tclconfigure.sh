# Auto-configure options for Mac OS X Universal build

# Written by Duncan Mortimer

macosx_universal_opts="-arch i386 -arch x86_64"
cflags="${cflags} ${macosx_universal_opts} -mmacosx-version-min=10.6"
cxxflags="${cxxflags} ${macosx_universal_opts} -mmacosx-version-min=10.6"
ldflags="${ldflags} -Wl,-search_paths_first ${macosx_universal_opts}"
configure_opts="${configure_opts} --disable-dependency-tracking"

