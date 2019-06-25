workflow "Push" {
  on = "push"
  resolves = ["Mirror to repo"]
}

workflow "Create ref" {
  on = "create"
  resolves = ["Possible dupe"]
}

workflow "Delete ref" {
  on = "delete"
  resolves = ["Possible dupe"]
}

workflow "Wiki-update" {
  on = "gollum"
  resolves = ["Mirror wiki"]
}

workflow "Nightly-mirror" {
  on = "schedule(49 7 * * *)" # 07:49 UTC daily
  resolves = ["Mirror to repo", "Mirror wiki"]
}

action "Mirror to repo" {
  uses = "docker://buildpack-deps:testing-scm"
  runs = ".github/main.workflow.sh"
  env = {
      MIRROR_URL = "git@github.com:nrc-fuels/MORFEUS-mirror.git"
  }
  secrets = ["IBB_PWLESS_DEPLOY_KEY", "GITHUB_TOKEN"]
}

action "Possible dupe" {
  needs = "Delay"
  uses = "docker://buildpack-deps:testing-scm"
  runs = ".github/main.workflow.sh"
  env = {
      MIRROR_URL = "git@github.com:nrc-fuels/MORFEUS-mirror.git"
  }
  secrets = ["IBB_PWLESS_DEPLOY_KEY", "GITHUB_TOKEN"]
}

action "Delay" {
  uses = "docker://alpine:latest"
  runs = ["sh", "-c", "sleep 30"]
}

action "Mirror wiki" {
  uses = "docker://buildpack-deps:testing-scm"
  runs = ".github/wiki-update.workflow.sh"
  env = {
      SOURCE_WIKI = "git@github.com:sourceryinstitute/MORFEUS-Source.wiki.git"
      MIRROR_WIKI = "git@github.com:nrc-fuels/MORFEUS-mirror.wiki.git"
  }
  secrets = ["IBB_PWLESS_DEPLOY_KEY", "GITHUB_TOKEN"]
}
