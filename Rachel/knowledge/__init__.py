"""Versioned, non-executable knowledge resources used by Rachel."""

from .profile import (
    KnowledgeProfile,
    KnowledgeProfileError,
    get_base_profile,
    resolve_pinned_knowledge_profile,
    resolve_knowledge_profile,
)

__all__ = [
    "KnowledgeProfile",
    "KnowledgeProfileError",
    "get_base_profile",
    "resolve_pinned_knowledge_profile",
    "resolve_knowledge_profile",
]
