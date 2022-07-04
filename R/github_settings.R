library(usethis)
use_git_config(user.name = "YOUR_GITHUB_USER",
               user.email = "YOUR_GITHUB_EMAIL")
git_vaccinate()
create_github_token()
gitcreds::gitcreds_set()
usethis::git_sitrep()
