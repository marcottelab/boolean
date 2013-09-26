require "mkmf"

$srcs = [
         'hypergeometric.c'
        ]


#if have_header("gsl/gsl_sf_exp.h", ["/usr/local/Cellar/gsl/1.15/include/"])
#  have_library("gsl")
#end


# Order matters here: ATLAS has to go after LAPACK: http://mail.scipy.org/pipermail/scipy-user/2007-January/010717.html
#$libs += " -llapack -lcblas -latlas "

$objs = %w{hypergeometric}.map { |i| i + ".o" }

CONFIG['CC'] = 'gcc'
#CONFIG['CXX'] = 'g++-4.7'

#dir_config("gsl", ["/usr/local/lib"])
find_library("gsl", "gsl_sf_lnfact")#, ["/usr/local/include/gsl/gsl_sf_gamma.h"])
find_header("gsl/gsl_sf_gamma.h", ["/usr/local/include/"])
# For release, these next two should both be changed to -O3.
$CFLAGS += " -O3 " #" -O0 -g "
# $CFLAGS += " -static -O0 -g "
$CPPFLAGS += " -O3 " #"-std=#{$CPP_STANDARD} " #" -O0 -g -std=#{$CPP_STANDARD} " #-fmax-errors=10 -save-temps
# $CPPFLAGS += " -static -O0 -g -std=#{$CPP_STANDARD} "

$CFLAGS += ' -std=c99 '
$libs   += ' -lgsl '

CONFIG['configure_args'].gsub!('-Wno-error=shorten-64-to-32', '')
CONFIG['CFLAGS'].gsub!('-Wno-error=shorten-64-to-32', '')
CONFIG['warnflags'].gsub!('-Wshorten-64-to-32', '') # doesn't work except in Mac-patched gcc (4.2)

#CONFIG['warnflags'].gsub!('-Wdeclaration-after-statement', '')
#CONFIG['warnflags'].gsub!('-Wimplicit-function-declaration', '')
# create_conf_h("hypergeometric_config.h")
create_makefile("hypergeometric")

