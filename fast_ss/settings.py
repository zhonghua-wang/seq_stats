from pydantic import BaseSettings


class Settings(BaseSettings):
    media_path: str = 'media'

    class Config:
        env_file = '.env'
        env_file_encoding = 'utf-8'
