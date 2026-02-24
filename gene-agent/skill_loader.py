"""
SkillLoader — reads .md skill files and parses YAML frontmatter.
Separated from agent.py so it can be imported independently.
"""

import os


class SkillLoader:
    def __init__(self, skills_dir: str = "./skills"):
        self.skills_dir = skills_dir
        self._cache = {}

    def load(self, skill_name: str) -> dict:
        if skill_name in self._cache:
            return self._cache[skill_name]

        skill_path = os.path.join(self.skills_dir, f"{skill_name}.md")
        if not os.path.exists(skill_path):
            raise FileNotFoundError(f"Skill file not found: {skill_path}")

        with open(skill_path, "r", encoding="utf-8") as f:
            raw = f.read()

        skill = self._parse(skill_name, raw)
        self._cache[skill_name] = skill
        return skill

    def _parse(self, name: str, raw: str) -> dict:
        metadata = {"name": name, "description": "", "version": "1.0", "tools": []}
        body = raw

        if raw.startswith("---"):
            parts = raw.split("---", 2)
            if len(parts) >= 3:
                frontmatter = parts[1].strip()
                body = parts[2].strip()
                for line in frontmatter.splitlines():
                    if ":" in line:
                        key, _, val = line.partition(":")
                        key = key.strip()
                        val = val.strip()
                        if key == "description":
                            metadata["description"] = val
                        elif key == "version":
                            metadata["version"] = val
                    elif line.strip().startswith("-"):
                        tool = line.strip().lstrip("-").strip()
                        if tool:
                            metadata["tools"].append(tool)

        return {"metadata": metadata, "body": body, "system_prompt": body}

    def list_skills(self) -> list:
        if not os.path.exists(self.skills_dir):
            return []
        return [
            f.replace(".md", "")
            for f in os.listdir(self.skills_dir)
            if f.endswith(".md")
        ]
