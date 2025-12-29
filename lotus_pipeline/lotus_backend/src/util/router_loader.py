import pkgutil
import importlib
import inspect
from fastapi import FastAPI, APIRouter


def auto_include_routers(app: FastAPI, package_name: str) -> None:
    """
    Automatically discover and register all APIRouter instances
    under the given package (Spring-style component scanning).
    """

    try:
        package = importlib.import_module(package_name)
    except ModuleNotFoundError as exc:
        print(f"Package not found: {package_name} ({exc})")
        return

    package_path = getattr(package, "__path__", None)
    if not package_path:
        return

    for _, module_name, _ in pkgutil.iter_modules(package_path):
        full_module_name = f"{package_name}.{module_name}"

        try:
            module = importlib.import_module(full_module_name)

            for _, obj in inspect.getmembers(module):
                if isinstance(obj, APIRouter):
                    print(f"Auto-register router: {full_module_name}")
                    app.include_router(obj)

        except Exception as exc:
            print(f"Failed to scan module {full_module_name}: {exc}")
