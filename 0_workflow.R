library(devtools)

check() # then...
# commit and push to github



# write a function, use roxygen skeleton to document it
# COMMON ROXYGEN TAGS
# @description @examples @examplesIf
# @export @family @inheritParams @param
# @returns @seealso @rdname
# ... then...
document()
use_package_doc()
check()

# commit and push to GH again


##### TRY INSTALLING FROM GITHUB
library(devtools)
devtools::install_github(repo = "https://github.com/joshua-nugent/tmle4rcts")
library(tmle4rcts)
?test_me_out

################ OR...
# skip it and just install freshest version
install()
