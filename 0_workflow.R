library(devtools)

check() # then...
# commit and push to github



# write a function, use roxygen skeleton to document it
# COMMON ROXYGEN TAGS
# @description @examples @examplesIf
# @export @family @inheritParams @param
# @returns @seealso @rdname
# To wrap blocks to be less then 80 characters wide: 
#  use Ctrl/Cmd + Shift + /
#  Or menu: code -> re-flow comment

# To other documentation:
#   \code{\link{function}}: function in this package.
# \code{\link[MASS]{abbey}}: function in another package.
# \link[=dest]{name}: link to dest, but show name.
# \code{\link[MASS:abbey]{name}}: link to function in another package, but show name.
# \linkS4class{abc}: link to an S4 class.
# To the web:
#   \url{http://rstudio.com}: a url.
# \href{http://rstudio.com}{Rstudio}:, a url with custom link text.
# \email{hadley@@rstudio.com} (note the doubled @): an email address.

# formatting with LATEX-style...
# \emph{will give italics}
# \strong{will give bold}
# \code{inline code format}
# \preformatted{for multiline code}
# \itemize{
    # \item First thing
    # \item Second thing
# }
# \tabular, etc.

# Special characters:
#   Use @@ to get @
#   Use \% to get %
#   Use \\ to get \

# use @inheritParams AnotherFunctionWithASharedArgument
# to copy parameters from that other function (not recursively tho)

# ... then...
devtools::document()
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
