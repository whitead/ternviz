from modal import Image, Volume, Stub, asgi_app, Secret
from fastapi import FastAPI, status, UploadFile, HTTPException, Depends
from fastapi.responses import Response
from fastapi.security import HTTPBearer, HTTPAuthorizationCredentials
import os
import uuid

VMD_URL = "https://storage.googleapis.com/fh-modal-artifacts/vmd-1.9.4a55.bin.LINUXAMD64-CUDA102-OptiX650-OSPRay185-RTXRTRT.opengl.tar.gz"
stub = Stub("ternviz-pdb")
vol = Volume.from_name("ternviz-scratch", create_if_missing=True)
volume_map = {"/var/tmp": vol}
web_app = FastAPI()
auth_scheme = HTTPBearer()

image = (
    Image.debian_slim(python_version="3.11")
    .apt_install("git", "wget", "libgomp1", "ffmpeg")
    .pip_install("git+https://github.com/whitead/ternviz.git")
    .run_commands(f"wget {VMD_URL} -O /tmp/vmd.tar.gz && tar -xvf /tmp/vmd.tar.gz -C /tmp")
    .run_commands("cd /tmp/vmd-* && ./configure && cd src && make install && hash")
)

with image.imports():
    import ternviz
    import os



@stub.function(image=image, cpu=1)
def get_pdb(pdb_query: str):
    pdb_id, p = ternviz.get_pdb(pdb_query)
    p.seek(0) # reset the file
    return pdb_id, p.read()


@stub.function(cpu=8, image=image, volumes=volume_map)
def pdb_render(pdb_id: str, pdb_file: str, start: int, stop: int):
    with open("/tmp/pdb.pdb", "wb") as f:
        f.write(pdb_file)
    ternviz.render(
        "/tmp/pdb.pdb",
        800,
        id=pdb_id,
        vmd="/usr/local/bin/vmd",
        script_name="render-pdb.vmd",
        color="black",
        args=[str(start), str(stop), "1.0", "Structure"],
    )
    vol.commit()
    return True

@stub.function(cpu=8, image=image, volumes=volume_map)
def movie(name: str, pdb_id: str):
    vol.reload()
    m = ternviz.movie(
            pdb_id,
            ffmpeg="ffmpeg",
            short_name=name,
            color="white",
    )
    # return the binary file
    with open(m, "rb") as f:
        result = f.read()
    for f in os.listdir("/var/tmp"):
        os.remove(f"/var/tmp/{f}")
    vol.commit()
    return result

@web_app.post("/render/")
async def render(pdb_file: UploadFile, message: str = "",  token: HTTPAuthorizationCredentials = Depends(auth_scheme)):
    if token.credentials != os.environ["AUTH_TOKEN"]:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Incorrect bearer token",
            headers={"WWW-Authenticate": "Bearer"},
        )
    
    frames = 360
    job_id = str(uuid.uuid4()).split("-")[0]

    contents = await pdb_file.read()

    def _render_gen():
        for i in range(0, frames, 20):
            yield job_id, contents, i, i + 20
    
    # render map takes a generator
    _ = list(pdb_render.starmap(_render_gen()))

    blob = movie.remote(message, job_id)
    return Response(content=blob, media_type="media/mp4", headers={"Content-Disposition": f"attachment; filename=ternviz.mp4"})

    
@stub.function(secrets=[Secret.from_name("web-auth-token")])
@asgi_app(custom_domains=["viz.proteincrow.ai"])
def app():
    return web_app