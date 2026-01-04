import uvicorn
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from contextlib import asynccontextmanager

from infrastructure.database.connection import get_default_db_manager
from infrastructure.database.database_setup import db_setup
from infrastructure.filesystem.logging import setup_logging
from infrastructure.workspace_context import workspace_path_manager
from util.router_loader import auto_include_routers


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application startup and shutdown lifecycle."""

    print("System starting up...")
    # Initialize workspace context
    workspace_path_manager.initialize()
    # Initialize logging
    setup_logging(debug_mode=True)

    # Initialize database connections
    db_setup(ensure_schema=True)

    yield

    print("System shutting down...")
    get_default_db_manager().close_connection()


app = FastAPI(
    title="Lotus Biological Agent API",
    description="Backend engine for single-cell analysis",
    version="0.1.0",
    lifespan=lifespan,
)

# Allow all origins for local desktop / WebView usage
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Routers will be registered here automatically
auto_include_routers(app, "controller")


@app.get("/health")
def health_check():
    """Health check endpoint for frontend readiness detection."""
    return {"status": "ok", "version": "0.1.0"}


if __name__ == "__main__":
    # For local development only.
    # In production, the server is started by PyWebView or an external launcher.
    uvicorn.run(
        "main:app",
        host="127.0.0.1",
        port=8888,
        reload=True,
    )
