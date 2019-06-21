workflow "Push" {
  on = "push"
  resolves = ["Mirror to repo"]
}

workflow "Create-ref" {
  on = "create"
  resolves = ["Mirror to repo"]
}

workflow "Delete-ref" {
  on = "delete"
  resolves = ["Mirror to repo"]
}

action "Mirror to repo" {
  uses = "docker://buildpack-deps:testing-scm"
  runs = ".github/main.workflow.sh"
  env = {
      MIRROR_URL = "git@github.com:nrc-fuels/MORFEUS-mirror.git"
  }
  secrets = ["MIRROR_DEPLOYMENT_KEY", "SI_BOT_KEY", "GITHUB_TOKEN"]
}
