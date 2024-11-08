# Installing myTAI on debian

Build the image:
`docker build -t mytai-debian-check .`

This may take a while, on my machine it took around 25 minutes to build the image.

Run the container:
`docker run -v $(pwd)/output:/output mytai-debian-check`

This will create a myTAI.Rcheck directory with the output from the `R CMD check` command inside of `/output`
