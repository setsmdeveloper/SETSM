#!/bin/sh

git_description=$(cd "$1" && git describe --always --tags --dirty)

echo "#define GIT_DESCRIPTION \"${git_description}\"" > "${2}/_git_description.h"
    if ! cmp --silent -- "${2}/git_description.h" "${2}/_git_description.h"; then
        cp "${2}/_git_description.h" "${2}/git_description.h";
    fi
rm "${2}/_git_description.h"
