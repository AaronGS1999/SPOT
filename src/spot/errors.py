class SpotError(Exception):
    pass


class ConfigError(SpotError):
    pass


class InputError(SpotError):
    pass


class ToolError(SpotError):
    pass
